///////////////////////////////////////////////////////////////////////////////
//
// File: TurbulenceModel.cpp
//
// ESQUELETO DE TESIS - no forma parte de Nektar++ upstream.
///////////////////////////////////////////////////////////////////////////////

#include "TurbulenceModel.h"

#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <cmath>

namespace Nektar
{

TurbulenceModelFactory &GetTurbulenceModelFactory()
{
    static TurbulenceModelFactory instance;
    return instance;
}

TurbulenceModel::TurbulenceModel(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const int spacedim)
    : m_session(pSession), m_fields(pFields), m_spacedim(spacedim)
{
    // VariableConverter es el objeto que sabe pasar de variables conservativas
    // a primitivas. Se construye igual que en ArtificialDiffusion.cpp:54.
    m_varConv = MemoryManager<VariableConverter>::AllocateSharedPtr(m_session,
                                                                    spacedim);

    m_session->LoadParameter("PrT", m_PrT, 0.9);

    SetUpFilterWidth();
}

/**
 * @brief Calcula el ancho de filtro Delta en cada punto de cuadratura.
 *
 * Implementación sugerida (opción (a) del comentario en el .h):
 *     Delta_e = h_e / (P_e + 1)
 * uniforme dentro de cada elemento.
 *
 * Cómo obtener h_e: el camino más robusto en malla no estructurada es a partir
 * del volumen del elemento,
 *     h_e = V_e^(1/dim)
 * y V_e se obtiene integrando 1 sobre el elemento:
 *     V_e = (*exp)[e]->Integral(Array<OneD,NekDouble>(nq_e, 1.0));
 *
 * Ventaja: funciona igual en tets, hexas, prismas y pirámides, y no depende de
 * que la malla sea isótropa. Desventaja: en elementos muy estirados (capa
 * límite) subestima la dirección larga. Si hacés casos con pared, mirá las
 * definiciones anisótropas de Delta (Deardorff, Scotti).
 */
void TurbulenceModel::SetUpFilterWidth()
{
    const int nPts = m_fields[0]->GetNpoints();
    m_deltaFilter  = Array<OneD, NekDouble>(nPts, 0.0);

    // TODO: implementar. Esqueleto del bucle:
    //
    //   for (int e = 0; e < m_fields[0]->GetExpSize(); ++e)
    //   {
    //       auto exp     = m_fields[0]->GetExp(e);
    //       int  offset  = m_fields[0]->GetPhys_Offset(e);
    //       int  nqe     = exp->GetTotPoints();
    //       int  P       = exp->GetBasisNumModes(0) - 1;
    //
    //       Array<OneD, NekDouble> ones(nqe, 1.0);
    //       NekDouble volume = exp->Integral(ones);
    //       NekDouble h      = std::pow(volume, 1.0 / m_spacedim);
    //       NekDouble delta  = h / (P + 1.0);
    //
    //       Vmath::Fill(nqe, delta, &m_deltaFilter[offset], 1);
    //   }
    //
    // NOTA: esto se calcula UNA vez en el constructor porque la malla no se
    // mueve. Si algún día usás malla móvil (ALE), hay que recalcularlo.
}

/**
 * @brief Tensor de deformación desviador y su norma.
 *
 *     S_ij   = 0.5 * ( du_i/dx_j + du_j/dx_i )
 *     S^d_ij = S_ij - (1/3) * delta_ij * S_kk
 *     |S|    = sqrt( 2 * S^d_ij * S^d_ij )
 *
 * EL PUNTO CRÍTICO ES LA RESTA DE LA TRAZA. Sin ella, en una onda de
 * compresión pura (sin corte, sin turbulencia) el modelo ve |S| grande y
 * genera mu_t. Ese es el error número uno al portar un modelo incompresible.
 *
 * Convención de índices asumida:
 *     derivatives[j][i][q] = d(u_i)/d(x_j) en el punto de cuadratura q
 * VERIFICALA contra NavierStokesCFE::GetViscousFluxVector antes de confiar.
 * Si la invertís, el tensor simétrico sale igual (por suerte) pero cualquier
 * modelo que use la parte antisimétrica (WALE) sale mal.
 *
 * Almacenamiento de Sij: solo las componentes independientes, en el orden
 *     2D: [S00, S01, S11]                          (3 componentes)
 *     3D: [S00, S01, S02, S11, S12, S22]           (6 componentes)
 */
void TurbulenceModel::GetDeviatoricStrainRate(
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &derivatives,
    Array<OneD, Array<OneD, NekDouble>> &Sij, Array<OneD, NekDouble> &SNorm)
{
    const int nPts = m_fields[0]->GetNpoints();

    // TODO: implementar.
    //
    //   1. Calcular la traza (dilatación):
    //        div = du/dx + dv/dy + dw/dz
    //
    //   2. Para cada par (i,j):
    //        S_ij = 0.5*(derivatives[j][i][q] + derivatives[i][j][q])
    //        si i == j:  S_ij -= div/3
    //
    //   3. SNorm[q] = sqrt( 2 * sum_ij S_ij^2 )
    //      OJO: en la suma, los términos fuera de la diagonal cuentan DOS
    //      veces porque solo guardaste el triángulo superior. Ese factor 2
    //      olvidado es un bug clásico que te da un mu_t sistemáticamente
    //      chico y difícil de detectar (el caso "casi anda").
    //
    // Escribí un test unitario con un campo de velocidad analítico donde
    // conozcas |S| a mano (p.ej. corte puro u = (y, 0, 0) => |S| = 1) ANTES
    // de conectarlo al solver.

    boost::ignore_unused(derivatives, Sij);
    Vmath::Fill(nPts, 0.0, SNorm, 1);
}

/**
 * @brief kappa_t = cp * mu_t / Pr_t
 */
void TurbulenceModel::v_GetTurbulentThermalConductivity(
    const Array<OneD, NekDouble> &muT, Array<OneD, NekDouble> &kappaT)
{
    const int nPts = muT.size();

    NekDouble gamma, gasConstant;
    m_session->LoadParameter("Gamma", gamma, 1.4);
    m_session->LoadParameter("GasConstant", gasConstant, 287.058);

    const NekDouble cp = gamma * gasConstant / (gamma - 1.0);

    // kappaT = (cp / PrT) * muT
    Vmath::Smul(nPts, cp / m_PrT, muT, 1, kappaT, 1);
}

/**
 * @brief Presión de submalla isotrópica. Cero por defecto.
 *
 * Si la implementás (modelo de Yoshizawa):
 *     p_sgs = (2/3) * rho * C_I * Delta^2 * |S|^2,   C_I ~ 0.0066
 * y hay que sumarla a la presión termodinámica en el tensor de tensiones.
 * Relevante solo a Mach turbulento alto; despreciable en la mayoría de los
 * casos subsónicos. Si la despreciás, DECILO en el escrito con la
 * justificación (Mach turbulento del caso), no la omitas en silencio.
 */
void TurbulenceModel::v_GetSGSPressure(
    [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &physfield,
    [[maybe_unused]] const Array<OneD, Array<OneD, Array<OneD, NekDouble>>>
        &derivatives,
    Array<OneD, NekDouble> &pSGS)
{
    Vmath::Fill(pSGS.size(), 0.0, pSGS, 1);
}

} // namespace Nektar
