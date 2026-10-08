///////////////////////////////////////////////////////////////////////////////
//
// File: SVVOperator.h
//
// ESQUELETO DE TESIS - no forma parte de Nektar++ upstream.
//
// Description: Clase base abstracta para operadores de Spectral Vanishing
//              Viscosity (SVV) en el solver compresible.
//
// ---------------------------------------------------------------------------
// CÓMO LEER ESTE ARCHIVO
//
// Está comentado con dos niveles:
//   //  ...        -> comentario normal, sobre el problema (SVV, física, diseño)
//   // [C++ N]     -> nota didáctica sobre la construcción de C++ en esa línea
//
// El modelo a imitar es:
//   solvers/CompressibleFlowSolver/ArtificialDiffusion/ArtificialDiffusion.h
// Abrilo en paralelo: este archivo es deliberadamente isomorfo a aquel.
///////////////////////////////////////////////////////////////////////////////

// [C++ 1] HEADER GUARD. Un .h puede ser incluido por muchos .cpp, y a veces
// dos veces en el mismo .cpp por caminos distintos. Sin esto, el compilador
// vería la clase declarada dos veces y fallaría. La convención de Nektar++ es
// NEKTAR_<CARPETA>_<ARCHIVO>. El nombre solo tiene que ser único en todo el
// proyecto; su valor no importa.
#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_SVVOPERATOR_BASE
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_SVVOPERATOR_BASE

#include <string>

// [C++ 2] #include con <> busca en los directorios de include configurados por
// CMake; con "" busca primero al lado de este archivo. Nektar++ usa <> para
// cosas de la librería y "" para archivos hermanos.
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.h>
#include <MultiRegions/ExpList.h>

// [C++ 3] Todo Nektar++ vive dentro de este namespace. Sirve para que dos
// librerías puedan tener una clase con el mismo nombre sin chocar.
namespace Nektar
{

// [C++ 4] DECLARACIÓN ADELANTADA (forward declaration). Le dice al compilador
// "existe una clase llamada SVVOperator, ya te voy a contar cómo es". Alcanza
// para poder declarar punteros a ella en las líneas siguientes.
class SVVOperator;

// [C++ 5] typedef = alias de tipo. shared_ptr es un puntero con conteo de
// referencias: cuando la última copia del puntero desaparece, el objeto se
// destruye solo. En Nektar++ prácticamente nunca vas a escribir 'new'/'delete'.
typedef std::shared_ptr<SVVOperator> SVVOperatorSharedPtr;

// [C++ 6] LA FACTORY. Este es el mecanismo central de extensibilidad de
// Nektar++ y conviene entenderlo bien porque lo vas a usar dos veces.
//
// NekFactory es un template. Los parámetros entre <> significan:
//   std::string    -> la CLAVE. Es el string que va a aparecer en el XML de
//                     sesión, p.ej. <SOLVERINFO><I PROPERTY="SVVType"
//                     VALUE="ExpKernel"/></SOLVERINFO>
//   SVVOperator    -> el TIPO BASE que la factory produce
//   los siguientes -> los ARGUMENTOS que recibe el constructor
//
// Efecto práctico: una clase derivada se "registra" con un string, y en tiempo
// de ejecución podés pedir GetSVVOperatorFactory().CreateInstance("ExpKernel", ...)
// sin que el solver sepa que esa clase existe. Agregar un kernel SVV nuevo NO
// requiere tocar el solver.
typedef LibUtilities::NekFactory<
    std::string, SVVOperator, const LibUtilities::SessionReaderSharedPtr &,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &, const int>
    SVVOperatorFactory;

// [C++ 7] Declaración de la función que devuelve el singleton de la factory.
// Se DECLARA acá y se DEFINE en el .cpp. Esa separación (declarar en .h,
// definir en .cpp) es el modelo de compilación de C++ y es lo primero que
// conviene tener claro.
SVVOperatorFactory &GetSVVOperatorFactory();

/**
 * @class SVVOperator
 * @brief Operador de Spectral Vanishing Viscosity para el solver compresible.
 *
 * CONTEXTO DEL PROBLEMA
 * ---------------------
 * SVV agrega disipación SOLO a los modos polinomiales altos de cada elemento,
 * dejando intactos los modos bajos. Formalmente se agrega un término
 *
 *     eps_SVV * d/dx ( Q(x) * du/dx )
 *
 * donde Q es un operador que en espacio modal multiplica el modo n por un
 * kernel Q_n que vale 0 para n < M (modos resueltos) y crece hacia 1 para
 * n -> P (modos más altos). M = m_cutoffRatio * P.
 *
 * Ventaja frente a la viscosidad artificial de captura de choque: preserva la
 * precisión espectral de los modos bajos, así que la disipación introducida
 * no contamina la solución resuelta. Es lo que lo hace apto como modelo
 * implícito de submalla (iLES).
 *
 * ESTADO EN NEKTAR++
 * ------------------
 * - El kernel por elemento YA EXISTE:
 *     StdRegions::StdExpansion::SVVLaplacianFilter(array, mkey)
 *     StdRegions::StdExpansion::ExponentialFilter(array, alpha, exponent, cutoff)
 *   con implementaciones concretas en StdQuadExp, StdHexExp, StdTriExp,
 *   StdTetExp, StdPrismExp, StdPyrExp, StdSegExp.
 * - El bucle sobre elementos YA EXISTE:
 *     MultiRegions::ExpList::ExponentialFilter (ExpList.cpp:2352)
 * - La orquestación existe SOLO para el solver incompresible, en
 *   VelocityCorrectionScheme::SetUpSVV() / AppendSVVFactors().
 * - Para el solver COMPRESIBLE no existe nada. Eso es lo que agrega esta clase.
 *
 * DIFERENCIA CLAVE CON EL CASO INCOMPRESIBLE (leer antes de codear)
 * -----------------------------------------------------------------
 * El solver incompresible es IMPLÍCITO: resuelve un Helmholtz por paso, y mete
 * SVV como un factor extra en la StdMatrixKey (ver AppendSVVFactors, donde
 * setea eFactorSVVCutoffRatio y eFactorSVVDiffCoeff). Es decir, SVV queda
 * embebido en la matriz.
 *
 * El solver compresible (CompressibleFlowSystem) es EXPLÍCITO en su forma más
 * usada: arma un RHS y avanza con Runge-Kutta. No hay matriz donde meter SVV.
 * Por lo tanto acá hay dos caminos posibles y hay que elegir uno a conciencia:
 *
 *   (A) SVV como término adicional en el RHS: calcular el laplaciano filtrado
 *       y sumarlo al residuo viscoso. Es lo consistente con la formulación
 *       original de SVV y lo que permite hablar de "viscosidad".
 *
 *   (B) SVV como filtro de post-paso: aplicar el kernel a los coeficientes
 *       después de cada paso temporal. Es mucho más simple y barato, pero
 *       técnicamente es un filtrado espectral, no SVV, y la disipación efectiva
 *       depende del paso de tiempo (lo cual hay que reportar honestamente).
 *
 * Este esqueleto está armado para (A), con DoSVV() operando sobre el RHS.
 * Si terminás eligiendo (B), la interfaz casi no cambia: cambia dónde se llama.
 */
class SVVOperator
{
public:
    // [C++ 8] IDIOMA "NVI" (Non-Virtual Interface). Fijate el patrón: el método
    // PÚBLICO no es virtual y lo único que hace es llamar al método PROTEGIDO
    // v_...() que sí es virtual. Ese es el estilo uniforme de todo Nektar++.
    //
    // ¿Para qué? Para poder agregar en el futuro chequeos, timers o
    // precondiciones en el método público sin tocar ninguna clase derivada.
    // Vos, cuando derivés, SOLO redefinís las funciones v_.
    //
    // [C++ 9] Definir el cuerpo acá adentro de la clase lo hace implícitamente
    // 'inline': el compilador puede reemplazar la llamada por su contenido.
    // Para un wrapper de una línea es lo correcto.

    /// Aplica SVV al residuo. inarray = campos conservativos, outarray = RHS
    /// al que se le SUMA la contribución SVV.
    void DoSVV(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
               Array<OneD, Array<OneD, NekDouble>> &outarray)
    {
        v_DoSVV(inarray, outarray);
    }

    /// Devuelve el campo de viscosidad SVV equivalente, en puntos de
    /// cuadratura. NO se usa para avanzar la solución: se usa para
    /// DIAGNÓSTICO. Es lo que te va a permitir, en supersónico, graficar por
    /// separado cuánta disipación viene de SVV, cuánta de la captura de choque
    /// y cuánta del modelo de turbulencia. Sin esto no vas a poder defender
    /// que el modelo hace lo que decís que hace.
    void GetSVVViscosity(const Array<OneD, Array<OneD, NekDouble>> &physfield,
                         Array<OneD, NekDouble> &muSVV)
    {
        v_GetSVVViscosity(physfield, muSVV);
    }

    /// Devuelve el kernel Q_n evaluado para un orden polinomial dado.
    /// Existe SOLO para poder testearlo unitariamente contra la fórmula
    /// analítica sin tener que correr un caso completo. Ver test nivel 1 en
    /// 00_PLAN_APRENDIZAJE_CPP.md.
    void GetKernel(const int nModes, Array<OneD, NekDouble> &kernel)
    {
        v_GetKernel(nModes, kernel);
    }

protected:
    // [C++ 10] MIEMBROS DE DATOS. La convención de Nektar++ es prefijo 'm_'.
    // Son 'protected' y no 'private' para que las clases derivadas puedan
    // usarlos directamente (así lo hace ArtificialDiffusion).

    /// Lector del archivo XML de sesión: de acá salen todos los parámetros.
    LibUtilities::SessionReaderSharedPtr m_session;

    /// Los campos del solver. Array<OneD, X> es el array con conteo de
    /// referencias de Nektar++ (ver LibUtilities/BasicUtils/SharedArray.hpp).
    /// OJO: copiar un Array<OneD,> NO copia los datos, comparte el buffer.
    /// Esto sorprende viniendo de numpy y es fuente clásica de bugs.
    Array<OneD, MultiRegions::ExpListSharedPtr> m_fields;

    /// Dimensión espacial.
    int m_spacedim;

    /// Fracción del orden polinomial por encima de la cual actúa el kernel.
    /// M = m_cutoffRatio * P. Valor habitual 0.75 (es el default del solver
    /// incompresible; ver VelocityCorrectionScheme.cpp:1110).
    NekDouble m_cutoffRatio;

    /// Magnitud de la viscosidad SVV.
    NekDouble m_diffCoeff;

    /// Si es true, se aplica SVV solo a las componentes de momento y no a
    /// densidad ni energía. Es una decisión de diseño NO trivial en
    /// compresible: filtrar densidad puede violar conservación de masa, y
    /// filtrar energía puede alterar el salto de choque. Dejalo configurable
    /// y documentá qué elegiste.
    bool m_applyToMomentumOnly;

    // [C++ 11] CONSTRUCTOR PROTEGIDO. Al no ser público, nadie puede escribir
    // 'SVVOperator x(...)' desde afuera. La única forma de crear uno es a
    // través de la factory. Eso es intencional.
    SVVOperator(const LibUtilities::SessionReaderSharedPtr &pSession,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                const int spacedim);

    // [C++ 12] DESTRUCTOR VIRTUAL. Regla de oro: si una clase tiene alguna
    // función virtual, su destructor DEBE ser virtual. Si no, destruir por
    // puntero a la base no llama al destructor de la derivada -> fuga de
    // memoria. '= default' le dice al compilador "generá el cuerpo obvio".
    virtual ~SVVOperator() = default;

    // [C++ 13] Función virtual CON implementación por defecto. Las derivadas
    // pueden redefinirla, pero no están obligadas. Acá tiene sentido: el
    // recorrido de elementos y la suma al residuo son iguales para todos los
    // kernels; lo único que cambia es el kernel en sí.
    virtual void v_DoSVV(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    virtual void v_GetSVVViscosity(
        const Array<OneD, Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &muSVV);

    // [C++ 14] FUNCIÓN VIRTUAL PURA: el '= 0' final. Significa "no hay
    // implementación acá, toda clase derivada ESTÁ OBLIGADA a proveerla".
    // La consecuencia es que SVVOperator es una CLASE ABSTRACTA: no se puede
    // instanciar. Es exactamente el rol que queremos: define el contrato,
    // no la política.
    //
    // Este es el único método que distingue un kernel SVV de otro.
    virtual void v_GetKernel(const int nModes,
                             Array<OneD, NekDouble> &kernel) = 0;
};

} // namespace Nektar

#endif
