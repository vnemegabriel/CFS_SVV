///////////////////////////////////////////////////////////////////////////////
//
// File: SVVOperator.cpp
//
// ESQUELETO DE TESIS - no forma parte de Nektar++ upstream.
//
// Description: Implementación de la clase base SVVOperator.
//
// Comparar con:
//   solvers/CompressibleFlowSolver/ArtificialDiffusion/ArtificialDiffusion.cpp
///////////////////////////////////////////////////////////////////////////////

// [C++ 15] El .cpp incluye su propio .h. Esto es una comprobación gratuita de
// que el .h es autosuficiente (compila solo, sin depender de que alguien haya
// incluido otra cosa antes).
#include "SVVOperator.h"

#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <StdRegions/StdMatrixKey.h>

namespace Nektar
{

/**
 * [C++ 16] SINGLETON DE MEYERS.
 *
 * 'static' dentro de una función significa: esta variable se construye UNA sola
 * vez, la primera vez que se ejecuta la función, y sobrevive hasta el final del
 * programa. Devolver una referencia (&) y no una copia es esencial: queremos
 * que todos vean LA MISMA factory, si no los registros se perderían.
 *
 * Es thread-safe desde C++11. Es el mismo patrón exacto que
 * GetArtificialDiffusionFactory() en ArtificialDiffusion.cpp:41.
 */
SVVOperatorFactory &GetSVVOperatorFactory()
{
    static SVVOperatorFactory instance;
    return instance;
}

/**
 * [C++ 17] CONSTRUCTOR CON LISTA DE INICIALIZACIÓN.
 *
 * Los ': m_session(pSession), ...' antes de la llave se llaman lista de
 * inicialización. Construyen los miembros DIRECTAMENTE con ese valor, en vez
 * de construirlos vacíos y después asignarles. Es más eficiente y, para
 * miembros const o referencias, es la única forma posible.
 *
 * OJO: los miembros se inicializan en el orden en que están DECLARADOS en el
 * .h, no en el orden que escribís acá. Poner un orden distinto es una fuente
 * clásica de warnings y de bugs sutiles. Mantené el mismo orden.
 */
SVVOperator::SVVOperator(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const int spacedim)
    : m_session(pSession), m_fields(pFields), m_spacedim(spacedim)
{
    // [C++ 18] LoadParameter tiene la firma
    //     LoadParameter(nombre, variable_destino, valor_por_defecto)
    // El segundo argumento se pasa por referencia NO const, así que la función
    // lo MODIFICA. Ese es el idioma "parámetro de salida", muy usado en este
    // codebase (viene de una época previa a los tipos de retorno baratos).
    //
    // Los defaults 0.75 y 0.1 son los del solver incompresible
    // (VelocityCorrectionScheme.cpp:1110 y :1069). Son un punto de partida
    // razonable, pero para el compresible los vas a tener que recalibrar:
    // documentá el barrido que hagas.
    m_session->LoadParameter("SVVCutoffRatio", m_cutoffRatio, 0.75);
    m_session->LoadParameter("SVVDiffCoeff", m_diffCoeff, 0.1);

    // DefinesSolverInfo/GetSolverInfo leen el bloque <SOLVERINFO> del XML.
    m_applyToMomentumOnly =
        m_session->DefinesSolverInfo("SVVMomentumOnly") &&
        m_session->GetSolverInfo("SVVMomentumOnly") == "True";

    // [C++ 19] ASSERTL0 es la macro de aserción de Nektar++ (nivel 0 = siempre
    // activa, incluso en Release). Usala para validar entrada del usuario.
    // Firma: ASSERTL0(condicion_que_debe_ser_verdadera, "mensaje si falla").
    ASSERTL0(m_cutoffRatio > 0.0 && m_cutoffRatio <= 1.0,
             "SVVCutoffRatio debe estar en (0, 1].");
    ASSERTL0(m_diffCoeff >= 0.0, "SVVDiffCoeff no puede ser negativo.");
}

/**
 * @brief Aplica SVV y SUMA el resultado a outarray.
 *
 * IMPORTANTE: la convención en este solver es que el operador SUMA su
 * contribución al residuo, no lo sobreescribe. Ver
 * ArtificialDiffusion::v_DoArtificialDiffusion (ArtificialDiffusion.cpp:69),
 * que termina con Vmath::Vadd sobre outarray. Respetá esa convención o vas a
 * borrar silenciosamente el término convectivo.
 */
void SVVOperator::v_DoSVV(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    // [C++ 20] '.size()' sobre un Array<OneD,> devuelve el número de elementos.
    // size_t es un entero sin signo; el codebase mezcla int y size_t según la
    // época del archivo. Elegí uno y sé consistente dentro de tu clase.
    const size_t nVariables = inarray.size();
    const size_t nCoeffs    = m_fields[0]->GetNcoeffs();
    const size_t nPts       = m_fields[0]->GetNpoints();

    // Índices de las variables de momento: 1..m_spacedim en la ordenación
    // conservativa de Nektar++ (rho, rho*u, rho*v, [rho*w], E).
    const size_t iFirst = m_applyToMomentumOnly ? 1 : 0;
    const size_t iLast =
        m_applyToMomentumOnly ? (1 + m_spacedim) : nVariables;

    // [C++ 21] Construcción de un array de trabajo. La sintaxis
    //   Array<OneD, NekDouble>(n, 0.0)
    // reserva n doubles y los inicializa a 0.0. Sin el segundo argumento la
    // memoria queda SIN inicializar (basura). Inicializá siempre salvo que
    // tengas una razón medida para no hacerlo.
    Array<OneD, NekDouble> tmpCoeffs(nCoeffs, 0.0);
    Array<OneD, NekDouble> tmpPhys(nPts, 0.0);

    for (size_t i = iFirst; i < iLast; ++i)
    {
        // ---------------------------------------------------------------
        // PASO 1: pasar de espacio físico (puntos de cuadratura) a espacio
        // modal (coeficientes). SVV es una operación intrínsecamente modal:
        // el kernel actúa sobre el número de modo.
        // ---------------------------------------------------------------
        m_fields[i]->FwdTrans(inarray[i], tmpCoeffs);

        // ---------------------------------------------------------------
        // PASO 2: aplicar el kernel elemento por elemento.
        //
        // ACÁ ESTÁ LA DECISIÓN DE DISEÑO PRINCIPAL. Tres opciones:
        //
        // (a) Reusar el filtro exponencial ya existente en la librería:
        //         m_fields[i]->ExponentialFilter(tmpPhys, alpha, expo, cutoff);
        //     (ver MultiRegions/ExpList.cpp:2352). Es la ruta más rápida a un
        //     resultado, pero opera en espacio físico y es filtrado, no SVV.
        //
        // (b) Reusar el laplaciano SVV de StdRegions:
        //         StdRegions::StdMatrixKey mkey(...);
        //         mkey.SetConstFactor(StdRegions::eFactorSVVCutoffRatio, m_cutoffRatio);
        //         mkey.SetConstFactor(StdRegions::eFactorSVVDiffCoeff, m_diffCoeff);
        //         (*exp)[e]->SVVLaplacianFilter(coeffs_del_elemento, mkey);
        //     Es lo que hace el solver incompresible vía AppendSVVFactors
        //     (VelocityCorrectionScheme.cpp:1221). Semánticamente correcto,
        //     requiere armar bien la StdMatrixKey.
        //
        // (c) Escribir tu propio bucle modal usando v_GetKernel(). Máximo
        //     control, necesario si tu tesis propone un kernel nuevo.
        //     Es lo que este esqueleto asume.
        //
        // Elegí y JUSTIFICÁ. El comité te lo va a preguntar.
        // ---------------------------------------------------------------

        // TODO: bucle sobre elementos. El patrón canónico (copiado de
        // ExpList::ExponentialFilter) es:
        //
        //   for (int e = 0; e < m_fields[i]->GetExpSize(); ++e)
        //   {
        //       int offset  = m_fields[i]->GetCoeff_Offset(e);
        //       int nModesE = m_fields[i]->GetExp(e)->GetNcoeffs();
        //
        //       Array<OneD, NekDouble> kernel(nModesE, 0.0);
        //       GetKernel(nModesE, kernel);   // <- polimorfismo: llama a la derivada
        //
        //       // multiplicar coeficiente a coeficiente por el kernel
        //       Vmath::Vmul(nModesE, kernel, 1,
        //                            tmpCoeffs + offset, 1,
        //                            tmpCoeffs + offset, 1);
        //   }
        //
        // [C++ 22] Ese 'tmpCoeffs + offset' es aritmética de punteros sobre el
        // Array: produce una VISTA que arranca en offset, sin copiar. El '1'
        // repetido es el "stride" (de a cuántos elementos avanzar); casi
        // siempre 1. Vmath::Vmul(n, x, incx, y, incy, z, incz) hace
        // z[k] = x[k]*y[k]. Toda la aritmética vectorial del codebase se
        // escribe así, no con bucles a mano: está vectorizada y probada.

        // ---------------------------------------------------------------
        // PASO 3: volver a espacio físico y SUMAR al residuo.
        // ---------------------------------------------------------------
        m_fields[i]->BwdTrans(tmpCoeffs, tmpPhys);

        // Escalar por el coeficiente SVV: tmpPhys = m_diffCoeff * tmpPhys
        Vmath::Smul(nPts, m_diffCoeff, tmpPhys, 1, tmpPhys, 1);

        // Sumar: outarray[i] = outarray[i] + tmpPhys
        Vmath::Vadd(nPts, outarray[i], 1, tmpPhys, 1, outarray[i], 1);
    }
}

/**
 * @brief Viscosidad SVV equivalente, para diagnóstico.
 *
 * Implementación por defecto: una estimación dimensional
 *     mu_SVV ~ m_diffCoeff * h / P
 * con h el tamaño local de elemento y P el orden. No es exacta —SVV no es una
 * viscosidad puntual— pero sirve para COMPARAR magnitudes contra mu_artificial
 * (captura de choque) y mu_t (turbulencia) en el mismo gráfico, que es lo que
 * necesitás para la parte supersónica.
 *
 * Si tu formulación permite una expresión mejor, redefinila en la derivada.
 */
void SVVOperator::v_GetSVVViscosity(
    [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &muSVV)
{
    // [C++ 23] [[maybe_unused]] silencia el warning de "parámetro no usado".
    // El parámetro está en la firma porque las derivadas SÍ lo van a necesitar.
    // Es la alternativa moderna a comentar el nombre del parámetro.

    // TODO: obtener h por elemento. Camino usual:
    //   m_fields[0]->GetExp(e)->GetGeom()->GetCoordim()
    // y el jacobiano vía GetMetricInfo()->GetJac(...). Mirar cómo lo hace
    // NonSmoothShockCapture::v_GetArtificialViscosity, que necesita lo mismo.
    Vmath::Fill(muSVV.size(), 0.0, muSVV, 1);
}

} // namespace Nektar
