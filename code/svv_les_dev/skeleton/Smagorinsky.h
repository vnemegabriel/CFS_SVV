///////////////////////////////////////////////////////////////////////////////
//
// File: Smagorinsky.h
//
// ESQUELETO DE TESIS - no forma parte de Nektar++ upstream.
//
// Description: Modelo de Smagorinsky de coeficiente constante, formulación
//              compresible (Favre). Es el modelo de referencia mínimo.
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_SMAGORINSKY
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_SMAGORINSKY

#include "TurbulenceModel.h"

namespace Nektar
{

/**
 * @brief Smagorinsky compresible.
 *
 *     mu_t = rho * (Cs * Delta)^2 * |S~^d|
 *
 * con |S~^d| = sqrt(2 * S~^d_ij * S~^d_ij) y S~^d el tensor de deformación
 * DESVIADOR de las velocidades filtradas de Favre.
 *
 * VALORES DE Cs: 0.1 - 0.2. El clásico teórico es 0.17 (Lilly, turbulencia
 * isótropa), pero en la práctica se usa 0.1 porque 0.17 sobre-disipa en flujos
 * con corte. En alto orden, además, Cs interactúa con la elección de Delta y
 * con la disipación numérica del esquema: NO copies un Cs de un código de
 * volúmenes finitos de segundo orden y esperes el mismo comportamiento.
 * Calibralo contra el Taylor-Green que ya tenés.
 *
 * LIMITACIONES QUE DEBÉS DECLARAR EN EL ESCRITO:
 *   - No se anula en flujo laminar: mu_t > 0 aunque no haya turbulencia, así
 *     que contamina la zona de transición.
 *   - No se anula en la pared (necesita Van Driest).
 *   - Es puramente disipativo: no puede representar backscatter.
 *   - Cerca de un choque, |S| es enorme y el modelo mete mu_t espuria.
 * Las tres primeras las resuelve Vreman o WALE. La cuarta necesita un
 * limitador explícito, y por eso este esqueleto lo incluye.
 */
class Smagorinsky : public TurbulenceModel
{
public:
    friend class MemoryManager<Smagorinsky>;

    static TurbulenceModelSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const int spacedim)
    {
        TurbulenceModelSharedPtr p =
            MemoryManager<Smagorinsky>::AllocateSharedPtr(pSession, pFields,
                                                          spacedim);
        return p;
    }

    static std::string className;

protected:
    void v_GetEddyViscosity(
        const Array<OneD, Array<OneD, NekDouble>> &physfield,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &derivatives,
        Array<OneD, NekDouble> &muT) override;

private:
    Smagorinsky(const LibUtilities::SessionReaderSharedPtr &pSession,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                const int spacedim);

    ~Smagorinsky() override = default;

    /// Constante de Smagorinsky.
    NekDouble m_Cs;

    /// Si es true, se apaga mu_t donde la dilatación es fuertemente negativa
    /// (compresión ~ choque). Ver punto 4 en TurbulenceModel.h.
    bool m_shockLimiter;

    /// Umbral de dilatación normalizada para el limitador.
    NekDouble m_dilatationThreshold;
};

} // namespace Nektar

#endif
