///////////////////////////////////////////////////////////////////////////////
//
// File: Smagorinsky.cpp
//
// ESQUELETO DE TESIS - no forma parte de Nektar++ upstream.
///////////////////////////////////////////////////////////////////////////////

#include "Smagorinsky.h"

#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <cmath>

namespace Nektar
{

std::string Smagorinsky::className =
    GetTurbulenceModelFactory().RegisterCreatorFunction(
        "Smagorinsky", Smagorinsky::create,
        "Modelo de Smagorinsky de coeficiente constante (compresible).");

Smagorinsky::Smagorinsky(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const int spacedim)
    : TurbulenceModel(pSession, pFields, spacedim)
{
    // 0.1 y no 0.17: ver la discusión en Smagorinsky.h. Calibralo.
    m_session->LoadParameter("Cs", m_Cs, 0.1);

    m_shockLimiter = m_session->DefinesSolverInfo("SGSShockLimiter") &&
                     m_session->GetSolverInfo("SGSShockLimiter") == "True";

    m_session->LoadParameter("SGSDilatationThreshold",
                             m_dilatationThreshold, 0.05);

    ASSERTL0(m_Cs >= 0.0, "Cs no puede ser negativo.");
}

/**
 * @brief mu_t = rho * (Cs*Delta)^2 * |S~^d|
 */
void Smagorinsky::v_GetEddyViscosity(
    const Array<OneD, Array<OneD, NekDouble>> &physfield,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &derivatives,
    Array<OneD, NekDouble> &muT)
{
    const int nPts = m_fields[0]->GetNpoints();

    // Número de componentes independientes del tensor simétrico.
    const int nSij = (m_spacedim == 3) ? 6 : 3;

    Array<OneD, Array<OneD, NekDouble>> Sij(nSij);
    for (int k = 0; k < nSij; ++k)
    {
        Sij[k] = Array<OneD, NekDouble>(nPts, 0.0);
    }
    Array<OneD, NekDouble> SNorm(nPts, 0.0);

    // Helper de la clase base: hace el trabajo pesado y la resta de la traza.
    GetDeviatoricStrainRate(derivatives, Sij, SNorm);

    // physfield[0] es la densidad en la ordenación conservativa de Nektar++.
    // VERIFICALO contra VariableConverter antes de confiar: si te equivocás de
    // índice, el código compila perfecto y da resultados sutilmente mal.
    const Array<OneD, const NekDouble> rho = physfield[0];

    // ------------------------------------------------------------------
    // mu_t = rho * (Cs*Delta)^2 * |S|
    //
    // Bucle explícito y no Vmath porque hay una operación por punto que
    // depende de m_deltaFilter[q]. Para claridad está escrito a mano; si el
    // profiling muestra que importa, se puede vectorizar precomputando
    // (Cs*Delta)^2 una sola vez en el constructor (que es lo correcto: Delta
    // no cambia en el tiempo).
    // ------------------------------------------------------------------
    for (int q = 0; q < nPts; ++q)
    {
        const NekDouble CsDelta = m_Cs * m_deltaFilter[q];
        muT[q] = rho[q] * CsDelta * CsDelta * SNorm[q];
    }

    // ------------------------------------------------------------------
    // LIMITADOR DE CHOQUE (ver punto 4 en TurbulenceModel.h).
    //
    // Idea: donde hay compresión fuerte (div(u) muy negativo comparado con
    // |S|), lo que el modelo "ve" es un choque, no turbulencia. Apagar mu_t
    // ahí evita ensanchar el choque artificialmente y evita contar dos veces
    // la disipación que ya mete ArtificialDiffusion.
    //
    // Criterio sugerido (tipo Ducros):
    //     theta = div(u)
    //     f = theta^2 / (theta^2 + |omega|^2 + eps)
    //     f -> 1 en choque, f -> 0 en turbulencia con vorticidad
    //     muT *= (1 - f)
    //
    // El sensor de Ducros es más robusto que un umbral duro sobre la
    // dilatación, y es el que se usa en la literatura de LES compresible con
    // choques. Vale la pena implementarlo bien: es un resultado defendible
    // por sí mismo en la tesis.
    // ------------------------------------------------------------------
    if (m_shockLimiter)
    {
        // TODO: implementar el sensor de Ducros usando derivatives.
        // Necesitás div(u) (traza, ya la calculaste en
        // GetDeviatoricStrainRate: considerá devolverla en vez de recalcular)
        // y |omega|^2 a partir de la parte antisimétrica del gradiente.
    }

    // Salvaguarda: mu_t no puede ser negativa. Con Smagorinsky puro nunca lo
    // es, pero si más adelante implementás el modelo dinámico SÍ puede serlo
    // (backscatter) y desestabiliza el solver. Dejá el clip acá desde ahora.
    for (int q = 0; q < nPts; ++q)
    {
        muT[q] = std::max(muT[q], 0.0);
    }
}

} // namespace Nektar
