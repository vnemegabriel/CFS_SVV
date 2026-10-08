///////////////////////////////////////////////////////////////////////////////
//
// File: TurbulenceModel.h
//
// ESQUELETO DE TESIS - no forma parte de Nektar++ upstream.
//
// Description: Clase base abstracta para modelos de turbulencia de submalla
//              (LES) en flujo compresible, válidos de subsónico a supersónico.
//
// A esta altura ya viste el patrón dos veces (ArtificialDiffusion y
// SVVOperator). Las notas [C++ N] acá son menos densas: el foco pasa a las
// decisiones de MODELADO, que es donde está la dificultad real de esta parte.
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_TURBULENCEMODEL_BASE
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_TURBULENCEMODEL_BASE

#include <string>

#include <CompressibleFlowSolver/Misc/VariableConverter.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>

namespace Nektar
{

class TurbulenceModel;

typedef std::shared_ptr<TurbulenceModel> TurbulenceModelSharedPtr;

typedef LibUtilities::NekFactory<
    std::string, TurbulenceModel, const LibUtilities::SessionReaderSharedPtr &,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &, const int>
    TurbulenceModelFactory;

TurbulenceModelFactory &GetTurbulenceModelFactory();

/**
 * @class TurbulenceModel
 * @brief Modelo de submalla para el CompressibleFlowSolver.
 *
 * ============================================================================
 * DÓNDE SE ENGANCHA ESTO (leer antes de escribir una línea)
 * ============================================================================
 *
 * En NavierStokesCFE la viscosidad molecular se calcula en:
 *
 *     NavierStokesCFE::GetViscosityAndThermalCondFromTemp(temperature,
 *                                                          mu, thermalCond)
 *     (declarada en NavierStokesCFE.h:155, con su versión kernel en :176)
 *
 * que a su vez llama a m_varConv->GetDynamicViscosity(temperature)
 * (NavierStokesCFE.h:192). Ese mu se usa después en GetViscousFluxVector.
 *
 * EL ENGANCHE LIMPIO es: después de calcular mu molecular, sumarle mu_t:
 *
 *     mu_efectiva  = mu_laminar + mu_t
 *     kappa_efect. = kappa_laminar + cp * mu_t / Pr_t
 *
 * Ventaja de este enganche: NO tocás el operador de difusión ni el flujo
 * viscoso. Todo el resto del solver sigue funcionando igual, solo que con una
 * viscosidad más grande donde el modelo lo pide. Es el cambio de menor
 * superficie posible, y por lo tanto el de menor riesgo de romper algo.
 *
 * ============================================================================
 * LO QUE HACE ESTO NO TRIVIAL EN COMPRESIBLE (sub -> supersónico)
 * ============================================================================
 *
 * 1. FILTRADO DE FAVRE. En compresible el filtrado es ponderado por densidad:
 *    u~ = <rho*u> / <rho>. Las fórmulas incompresibles de Smagorinsky escritas
 *    con rho constante están MAL en supersónico. En la práctica esto aparece
 *    como el factor rho que multiplica a nu_t:
 *        mu_t = rho * (Cs*Delta)^2 * |S~|
 *    Notá que devolvemos mu_t (dinámica), NO nu_t (cinemática). Con rho
 *    variable la distinción importa.
 *
 * 2. LA PARTE DILATACIONAL DEL TENSOR. En compresible div(u) != 0. El tensor
 *    de deformación que entra en el modelo debe ser el DESVIADOR:
 *        S~_ij^d = S~_ij - (1/3) * delta_ij * S~_kk
 *    Si usás S~_ij sin quitarle la traza, en una onda de compresión el modelo
 *    ve una "deformación" enorme que no es turbulencia y te mete viscosidad
 *    turbulenta espuria justo en el choque.
 *
 * 3. TENSIÓN ISOTRÓPICA DE SUBMALLA (término de Yoshizawa). Estrictamente
 *    tau_kk != 0 en compresible y aporta una presión de submalla
 *        p_sgs = (2/3) * rho * C_I * Delta^2 * |S~|^2
 *    Muchos trabajos la desprecian porque es pequeña a Mach turbulento bajo.
 *    A Mach turbulento alto deja de serlo. DECIDÍ EXPLÍCITAMENTE y justificá;
 *    el esqueleto deja el hook para calcularla.
 *
 * 4. INTERACCIÓN CON LA CAPTURA DE CHOQUE. El solver ya tiene
 *    ArtificialDiffusion metiendo mu_artificial cerca de los choques. Con LES
 *    activo hay dos (o tres, contando SVV) mecanismos disipativos superpuestos.
 *    Riesgos: doble contabilización de disipación, y que Smagorinsky trate el
 *    choque como capa de corte. Mitigaciones a evaluar:
 *      - usar un modelo con mejor comportamiento cerca de la pared y de
 *        discontinuidades (Vreman, WALE, sigma) en vez de Smagorinsky puro;
 *      - apagar mu_t donde el sensor de choque esté activo;
 *      - limitar por dilatación: si div(u) << 0 (compresión fuerte), reducir mu_t.
 *    El método GetEddyViscosity recibe el sensor por eso.
 *
 * 5. COMPORTAMIENTO EN LA PARED. Smagorinsky estándar da mu_t != 0 en la pared,
 *    donde debería anularse como y^3. Si vas a hacer casos con pared (que es
 *    donde vive la mayor parte de la aerodinámica interesante), Smagorinsky
 *    puro necesita amortiguamiento de Van Driest, o directamente usá WALE /
 *    Vreman que se anulan solos. En Taylor-Green (sin paredes) no importa;
 *    en un perfil alar sí.
 *
 * ============================================================================
 * ORDEN DE IMPLEMENTACIÓN SUGERIDO
 * ============================================================================
 *   1. Smagorinsky constante  -> es el "hola mundo", sirve para validar toda
 *                                la plomería con un caso Taylor-Green.
 *   2. Vreman o WALE          -> los que realmente vas a defender: se anulan
 *                                en flujo laminar y cerca de la pared.
 *   3. Smagorinsky dinámico   -> solo si el cronograma lo permite; requiere
 *                                filtro de test y promediado, que en mallas no
 *                                estructuradas de alto orden es un problema en
 *                                sí mismo. No lo prometas en la propuesta.
 */
class TurbulenceModel
{
public:
    /// Calcula la viscosidad turbulenta mu_t en los puntos de cuadratura.
    ///
    /// @param physfield   Variables físicas (rho, u, v, [w], p/T) en puntos
    ///                    de cuadratura.
    /// @param derivatives derivatives[j][i] = d(u_i)/d(x_j). El solver ya las
    ///                    calcula para el flujo viscoso: pasalas, no las
    ///                    recalcules (es la parte cara).
    /// @param muT         SALIDA: viscosidad turbulenta DINÁMICA (incluye rho).
    void GetEddyViscosity(
        const Array<OneD, Array<OneD, NekDouble>> &physfield,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &derivatives,
        Array<OneD, NekDouble> &muT)
    {
        v_GetEddyViscosity(physfield, derivatives, muT);
    }

    /// Presión de submalla isotrópica (término de Yoshizawa). Ver punto 3.
    /// La implementación por defecto devuelve cero, que es la hipótesis
    /// habitual a Mach turbulento bajo.
    void GetSGSPressure(
        const Array<OneD, Array<OneD, NekDouble>> &physfield,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &derivatives,
        Array<OneD, NekDouble> &pSGS)
    {
        v_GetSGSPressure(physfield, derivatives, pSGS);
    }

    /// Conductividad térmica turbulenta a partir de mu_t vía el número de
    /// Prandtl turbulento: kappa_t = cp * mu_t / Pr_t.
    void GetTurbulentThermalConductivity(const Array<OneD, NekDouble> &muT,
                                         Array<OneD, NekDouble> &kappaT)
    {
        v_GetTurbulentThermalConductivity(muT, kappaT);
    }

protected:
    LibUtilities::SessionReaderSharedPtr m_session;
    Array<OneD, MultiRegions::ExpListSharedPtr> m_fields;
    VariableConverterSharedPtr m_varConv;
    int m_spacedim;

    /// Prandtl turbulento. 0.9 es el valor convencional en aerodinámica; en
    /// capas de mezcla libres se usa 0.5-0.7. Hacelo configurable y reportá
    /// cuál usaste: es una de las primeras cosas que te van a preguntar.
    NekDouble m_PrT;

    /// Longitud de filtro Delta por punto de cuadratura. En alto orden la
    /// elección NO es obvia: las tres opciones habituales son
    ///   (a) Delta = h_elemento / (P+1)     <- la más usada, "espaciado modal"
    ///   (b) Delta = V_elemento^(1/3) / (P+1)
    ///   (c) Delta = distancia local entre puntos de cuadratura (varía dentro
    ///       del elemento porque los puntos GLL se agrupan en los bordes)
    /// Afecta a mu_t como Delta^2, o sea que la elección importa MUCHO.
    /// Documentá cuál usaste; idealmente mostrá sensibilidad a esta elección.
    Array<OneD, NekDouble> m_deltaFilter;

    TurbulenceModel(const LibUtilities::SessionReaderSharedPtr &pSession,
                    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                    const int spacedim);

    virtual ~TurbulenceModel() = default;

    /// Único método obligatorio: es lo que define al modelo.
    virtual void v_GetEddyViscosity(
        const Array<OneD, Array<OneD, NekDouble>> &physfield,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &derivatives,
        Array<OneD, NekDouble> &muT) = 0;

    virtual void v_GetSGSPressure(
        const Array<OneD, Array<OneD, NekDouble>> &physfield,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &derivatives,
        Array<OneD, NekDouble> &pSGS);

    virtual void v_GetTurbulentThermalConductivity(
        const Array<OneD, NekDouble> &muT, Array<OneD, NekDouble> &kappaT);

    /// Helper compartido: calcula el tensor de deformación DESVIADOR
    /// S~_ij - (1/3) delta_ij S~_kk y su norma |S~| = sqrt(2 S~_ij S~_ij).
    /// Está en la base porque TODOS los modelos algebraicos lo necesitan
    /// (Smagorinsky, Vreman, WALE, sigma). Escribirlo una vez, bien testeado,
    /// te ahorra el bug más común de esta parte de la tesis.
    void GetDeviatoricStrainRate(
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &derivatives,
        Array<OneD, Array<OneD, NekDouble>> &Sij,
        Array<OneD, NekDouble> &SNorm);

    /// Calcula m_deltaFilter. Llamado una vez desde el constructor.
    void SetUpFilterWidth();
};

} // namespace Nektar

#endif
