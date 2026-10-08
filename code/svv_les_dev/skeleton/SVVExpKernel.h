///////////////////////////////////////////////////////////////////////////////
//
// File: SVVExpKernel.h
//
// ESQUELETO DE TESIS - no forma parte de Nektar++ upstream.
//
// Description: Kernel SVV exponencial (Maday-Tadmor). Clase derivada concreta.
//
// Comparar con:
//   solvers/CompressibleFlowSolver/ArtificialDiffusion/NonSmoothShockCapture.h
//
// Este archivo es el MOLDE que vas a repetir para cada kernel nuevo. Si mañana
// querés un kernel de potencia o el kernel DG, copiás este archivo, cambiás
// tres líneas y v_GetKernel. Nada más.
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_SVVEXPKERNEL
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_SVVEXPKERNEL

#include "SVVOperator.h"

namespace Nektar
{

/**
 * @brief Kernel SVV exponencial.
 *
 * FÍSICA / NUMÉRICA
 * -----------------
 * Para el modo n de un elemento con P+1 modos, con corte M = ratio * P:
 *
 *     Q_n = 0                                    si n <= M
 *     Q_n = exp( -(n - P)^2 / (n - M)^2 )        si n > M
 *
 * Esta es la forma usada en Nektar++ para el solver incompresible (ver el
 * mensaje "SVV (... Exp Kernel(cut-off = ...)" en
 * VelocityCorrectionScheme.cpp:593). Propiedades relevantes:
 *   - Q_M = 0 exactamente: los modos resueltos no se tocan.
 *   - Q_P = 1: el modo más alto recibe disipación plena.
 *   - Transición suave: no introduce oscilaciones de Gibbs propias.
 *
 * Es el kernel por defecto razonable. Los otros dos que vale la pena
 * implementar después, ya que existen en la rama incompresible y te dan una
 * comparación gratis, son el "power kernel" y el "DG kernel".
 */

// [C++ 24] HERENCIA PÚBLICA. 'class Derivada : public Base' significa
// "toda SVVExpKernel ES-UN SVVOperator". Consecuencia práctica: podés guardar
// una SVVExpKernel en un SVVOperatorSharedPtr y el resto del solver nunca
// necesita saber cuál kernel es. Eso es polimorfismo, y es todo el punto del
// diseño.
class SVVExpKernel : public SVVOperator
{
public:
    // [C++ 25] 'friend' le da a MemoryManager acceso a nuestro constructor
    // privado. Sin esto, MemoryManager no podría construir la clase. Es
    // boilerplate: copialo tal cual, igual que en NonSmoothShockCapture.h:50.
    friend class MemoryManager<SVVExpKernel>;

    // [C++ 26] FUNCIÓN MIEMBRO STATIC. 'static' acá significa que NO necesita
    // un objeto para ser llamada: es SVVExpKernel::create(...), no
    // unObjeto.create(...). Tiene que ser static porque la factory la llama
    // justamente para crear el primer objeto — todavía no hay ninguno.
    //
    // Fijate que DEVUELVE un SVVOperatorSharedPtr (el tipo BASE), no un
    // puntero a SVVExpKernel. Eso es lo que permite que la factory tenga un
    // tipo de retorno uniforme para todos los kernels.
    static SVVOperatorSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const int spacedim)
    {
        // MemoryManager es el asignador de Nektar++. AllocateSharedPtr
        // construye el objeto y te devuelve un shared_ptr ya listo. Los
        // argumentos se reenvían al constructor.
        SVVOperatorSharedPtr p =
            MemoryManager<SVVExpKernel>::AllocateSharedPtr(pSession, pFields,
                                                           spacedim);
        return p;
    }

    // [C++ 27] VARIABLE MIEMBRO STATIC. Existe UNA sola para toda la clase, no
    // una por objeto. Acá se DECLARA; se DEFINE (y se le da valor) en el .cpp.
    // El truco es que su inicialización en el .cpp es lo que dispara el
    // registro en la factory antes de que arranque main().
    static std::string className;

protected:
    // [C++ 28] 'override' no es obligatorio pero SIEMPRE ponelo. Le pide al
    // compilador que verifique que efectivamente estás redefiniendo un virtual
    // de la base. Si te equivocás en un 'const' o en un tipo, en vez de crear
    // silenciosamente una función nueva que nunca se llama (bug clásico y
    // dificilísimo de encontrar), el compilador te frena.
    void v_GetKernel(const int nModes,
                     Array<OneD, NekDouble> &kernel) override;

private:
    // Constructor privado: solo create() y MemoryManager pueden invocarlo.
    SVVExpKernel(const LibUtilities::SessionReaderSharedPtr &pSession,
                 const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                 const int spacedim);

    ~SVVExpKernel() override = default;
};

} // namespace Nektar

#endif
