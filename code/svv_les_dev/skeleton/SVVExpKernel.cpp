///////////////////////////////////////////////////////////////////////////////
//
// File: SVVExpKernel.cpp
//
// ESQUELETO DE TESIS - no forma parte de Nektar++ upstream.
///////////////////////////////////////////////////////////////////////////////

#include "SVVExpKernel.h"

#include <cmath>

namespace Nektar
{

/**
 * [C++ 29] EL REGISTRO EN LA FACTORY. Esta línea es la más importante del
 * archivo y la más críptica. Leída de adentro hacia afuera:
 *
 *   GetSVVOperatorFactory()          -> obtiene el singleton
 *   .RegisterCreatorFunction(a, b, c) -> le agrega una entrada
 *       a = "ExpKernel"              -> la clave que va en el XML de sesión
 *       b = SVVExpKernel::create     -> puntero a la función que sabe crearlo
 *       c = "..."                    -> descripción humana
 *
 * ¿Por qué está escrito como la inicialización de una variable static? Porque
 * las variables static a nivel de namespace se inicializan ANTES de main().
 * Es un truco para ejecutar código de registro sin que nadie lo llame
 * explícitamente. Por eso agregar un kernel nuevo no requiere tocar ni una
 * línea del solver: basta con compilar el archivo.
 *
 * TRAMPA CONOCIDA: si linkeás esto en una librería estática y nadie referencia
 * el símbolo, el linker puede descartar el objeto y el registro nunca ocurre.
 * En Nektar++ los solvers se compilan de forma que esto funcione, pero si tu
 * clase "no aparece" en tiempo de ejecución sin ningún error, es casi siempre
 * esto.
 */
std::string SVVExpKernel::className =
    GetSVVOperatorFactory().RegisterCreatorFunction(
        "ExpKernel", SVVExpKernel::create,
        "Kernel SVV exponencial de Maday-Tadmor.");

/**
 * [C++ 30] El constructor de la derivada DEBE construir la base primero. Eso
 * es lo que hace ': SVVOperator(pSession, pFields, spacedim)'. Si no lo
 * escribís, C++ intenta llamar al constructor sin argumentos de la base — que
 * acá no existe — y falla la compilación.
 *
 * El cuerpo queda vacío porque toda la lectura de parámetros comunes ya la
 * hizo la base. Si este kernel tuviera parámetros propios, se leerían acá.
 */
SVVExpKernel::SVVExpKernel(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const int spacedim)
    : SVVOperator(pSession, pFields, spacedim)
{
}

/**
 * @brief Kernel exponencial.
 *
 *     Q_n = 0                                si n <= M
 *     Q_n = exp( -(n-P)^2 / (n-M)^2 )        si n > M
 *
 * con P = nModes - 1 (orden polinomial) y M = round(m_cutoffRatio * P).
 *
 * ESTA es la única función que tenés que escribir para tener un kernel nuevo.
 * Todo lo demás es infraestructura compartida.
 */
void SVVExpKernel::v_GetKernel(const int nModes, Array<OneD, NekDouble> &kernel)
{
    ASSERTL1(kernel.size() >= static_cast<size_t>(nModes),
             "El array kernel es mas chico que nModes.");

    // [C++ 31] 'const' en variables locales: no cambia el rendimiento, pero
    // documenta la intención y evita que las modifiques por accidente 40
    // líneas más abajo. Usalo siempre que puedas.
    const int P = nModes - 1;
    const int M = static_cast<int>(m_cutoffRatio * P);

    // [C++ 32] static_cast<T>(x) es la conversión explícita de C++. Preferila
    // al casteo estilo C '(int)x': es buscable con grep y el compilador
    // verifica que la conversión sea legítima.

    for (int n = 0; n < nModes; ++n)
    {
        if (n <= M)
        {
            // Modos resueltos: intactos. Este es EL punto de SVV.
            kernel[n] = 0.0;
        }
        else
        {
            // Caso borde: si M == P el denominador se anula. Puede pasar con
            // cutoffRatio = 1.0 o con P chico. Manejalo explícitamente en vez
            // de dejar que salga un NaN que después vas a perseguir por horas.
            if (n == M)
            {
                kernel[n] = 0.0;
                continue;
            }

            const NekDouble num = static_cast<NekDouble>((n - P) * (n - P));
            const NekDouble den = static_cast<NekDouble>((n - M) * (n - M));

            kernel[n] = std::exp(-num / den);
        }
    }

    // NOTA PARA EL TEST UNITARIO (nivel 1 del plan):
    //   - kernel[n] debe ser exactamente 0.0 para todo n <= M
    //   - kernel[P] debe ser exactamente 1.0  (porque n-P = 0 => exp(0) = 1)
    //   - la secuencia debe ser monótona creciente en n > M
    // Esas tres propiedades se testean sin correr ningún caso de CFD y te
    // cubren la mayoría de los errores de implementación.
}

} // namespace Nektar
