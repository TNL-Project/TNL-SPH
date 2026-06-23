"""
sph_codegen/gen.py  –  AUTO-GENERATED file emitter

Generates 4 C++ header files from the Python scheme definition:
  - Variables.h       (fluid/boundary variable containers)
  - control.h          (model parameter configuration class)
  - Interactions.h     (class declaration with all required types/methods)
  - Interactions.hpp   (CUDA kernel implementations)

The generated code matches the hand-written WCSPH-DBC structure so it
plugs directly into SPHMultiset_CFD< Model >.
"""
from __future__ import annotations
from pathlib import Path
from .dsl  import (_B, _U, _C, Fn, Dot, Ternary, Lit, ConstRef, FieldRef, LocalRef, Stmt, LocalDecl)
from .emit import emit
from .scheme import _FieldDesc

# ─── entry ──────────────────────────────────────────────────────────────────
def generate_all(fluid_cls, boundary_cls, constants_cls,
                 interactions_cls, integration_cls=None,
                 output_dir=".", namespace="WCSPH_DBC"):
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    files = {
        "Variables.h":       _gen_variables(fluid_cls, boundary_cls),
        "control.h":         _gen_control(constants_cls, namespace),
        "Interactions.h":    _gen_interactions_h(interactions_cls, namespace),
        "Interactions.hpp":  _gen_interactions_hpp(interactions_cls, fluid_cls, boundary_cls, constants_cls, namespace),
    }
    for fname, src in files.items():
        p = out / fname
        p.write_text(src)
        print(f"  ✓  {fname:<28} ({src.count(chr(10))} lines)  → {p}")
    return files

# ─── helpers ─────────────────────────────────────────────────────────────────
def _fields(cls): return [(k,v) for k,v in vars(cls).items() if isinstance(v,_FieldDesc)]
def _field_names(cls): return [k for k,v in vars(cls).items() if isinstance(v,_FieldDesc)]

def _ctype(fd):
    """C++ type for a field descriptor."""
    return "VectorArrayType" if fd.dtype == "vector" else "ScalarArrayType"

def _ctype_short(fd):
    """Short C++ type (for lambda parameters)."""
    return "VectorType" if fd.dtype == "vector" else "RealType"

def _init_val(fd):
    """Initial value for a field in the particle loop."""
    return "VectorType( 0.f )" if fd.dtype == "vector" else "0.f"

def _view_suffix(pset):
    """View name suffix for a particle set."""
    return "" if pset == "fluid" else "_bound"

def _search_token(owner, target):
    """Search token variable name."""
    return "searchInFluid" if target == "fluid" else "searchInBound"

def _output_fields(pair):
    """Unique list of output field names in an interaction record."""
    seen = []
    for s in pair.stmts:
        if s.field not in seen:
            seen.append(s.field)
    return seen

def _phase_map(ic):
    """Group interaction records by phase name."""
    g = {}
    for r in ic._interactions:
        pn = "interaction" if r.set_i == "fluid" else "updateSolidBoundary"
        g.setdefault(pn, []).append(r)
    return g

def _phase_names(ic): return list(_phase_map(ic).keys())

# ─── collect all field names accessed with 'j' in a record ──────────────────
def _j_fields(record, fluid_cls, boundary_cls):
    """Find all field names accessed with j index in the record's expressions."""
    found = set()
    def walk(node):
        if isinstance(node, FieldRef):
            if node.idx == "j":
                found.add(node.name)
        if hasattr(node, '__dict__'):
            for v in vars(node).values():
                if isinstance(v, list):
                    for item in v: walk(item)
                elif hasattr(v, '__dict__') or hasattr(v, 'name'):
                    walk(v)
        # dataclass fields
        if hasattr(node, 'left'): walk(node.left)
        if hasattr(node, 'right'): walk(node.right)
        if hasattr(node, 'operand'): walk(node.operand)
        if hasattr(node, 'cond'): walk(node.cond)
        if hasattr(node, 'yes'): walk(node.yes)
        if hasattr(node, 'no'): walk(node.no)
        if hasattr(node, 'a'): walk(node.a)
        if hasattr(node, 'b'): walk(node.b)
        if hasattr(node, 'args') and isinstance(node.args, list):
            for a in node.args: walk(a)
        if hasattr(node, 'expr'): walk(node.expr)
    for s in record.stmts:
        walk(s.expr)
    for d in record.locals:
        walk(d.expr)
    return found

# ─── Variables.h ─────────────────────────────────────────────────────────────
def _gen_variables(fluid_cls, boundary_cls) -> str:
    L = ["// Variables.h  –  AUTO-GENERATED",
         "#pragma once", "",
         '#include "../../SPHTraits.h"',
         '#include <TNL/Particles/details/thrustExecPolicySelector.h>',
         '#include <thrust/gather.h>',
         '#include "../WCSPH_DBC/BoundaryConditionsTypes.h"', "",
         "namespace TNL {", "namespace SPH {", ""]
    L += _var_class(fluid_cls, "FluidVariables")
    L.append("")
    L += _boundary_var_class(boundary_cls)
    L.append("")
    # OpenBoundaryVariables (needed by framework, inherits FluidVariables)
    L += ["template< typename SPHState >",
          "class OpenBoundaryVariables : public FluidVariables< SPHState >",
          "{",
          "public:",
          "   using BaseType = FluidVariables< SPHState >;",
          "   using SPHTraitsType = typename BaseType::SPHTraitsType;",
          "   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;",
          "   using IndexArrayType = typename SPHTraitsType::IndexArrayType;",
          "",
          "   void setSize( const GlobalIndexType& size )",
          "   {",
          "      BaseType::setSize( size );",
          "      particleMark.setSize( size );",
          "      receivingParticleMark.setSize( size );",
          "   }",
          "",
          "   typename SPHTraitsType::IndexArrayType particleMark;",
          "   typename SPHTraitsType::IndexArrayType receivingParticleMark;",
          "};",
          ""]
    L += ["} // SPH", "} // TNL", "", "// end Variables.h"]
    return "\n".join(L)

def _var_class(cls, tname):
    fds = _fields(cls)
    L  = [f"template< typename SPHState >",
          f"class {tname}", "{", "public:",
          "   using SPHConfig = typename SPHState::SPHConfig;",
          "   using SPHTraitsType = SPHFluidTraits< SPHConfig >;",
          "   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;",
          "   using RealType = typename SPHTraitsType::RealType;",
          "   using ScalarArrayType = typename SPHTraitsType::ScalarArrayType;",
          "   using VectorArrayType = typename SPHTraitsType::VectorArrayType;",
          "   using IndexArrayType = typename SPHTraitsType::IndexArrayType;",
          "   using IndexArrayTypePointer = typename Pointers::SharedPointer< IndexArrayType, typename SPHConfig::DeviceType >;",
          ""]
    # Field declarations
    for n, fd in fds:
        arr = _ctype(fd)
        L.append(f"   {arr} {n};")
        if fd.t >= 2:
            L.append(f"   {arr} {n}_swap;")
    L.append("")
    # setSize
    L += ["   void setSize( const GlobalIndexType& size )", "   {"]
    for n, fd in fds:
        L.append(f"      {n}.setSize( size );")
        if fd.t >= 2:
            L.append(f"      {n}_swap.setSize( size );")
    L += ["   }", ""]
    # sortVariables
    sortable = [(n, fd) for n, fd in fds if fd.t >= 2]
    if sortable:
        L += ["   template< typename ParticlesPointer >",
              "   void sortVariables( ParticlesPointer& particles )", "   {"]
        for n, _ in sortable:
            L.append(f"      particles->reorderArray( {n}, {n}_swap );")
        L += ["   }", ""]
    # readVariables
    readable = [(n, fd) for n, fd in fds if fd.read]
    if readable:
        L += ["   template< typename ReaderType >",
              "   void readVariables( ReaderType& reader )", "   {"]
        for n, fd in readable:
            cap = _vtk_name(n)
            if fd.dtype == "scalar":
                L.append(f'      reader.template readParticleVariable< ScalarArrayType, typename ScalarArrayType::ValueType >( {n}, "{cap}" );')
            else:
                L.append("      if constexpr( SPHConfig::spaceDimension == 2 )")
                L.append(f'         reader.template readParticleVariable2D< VectorArrayType, typename VectorArrayType::ValueType::ValueType >( {n}, "{cap}" );')
                L.append("      if constexpr( SPHConfig::spaceDimension == 3 )")
                L.append(f'         reader.template readParticleVariable3D< VectorArrayType, typename VectorArrayType::ValueType::ValueType >( {n}, "{cap}" );')
        L += ["   }", ""]
    # writeVariables
    writable = [(n, fd) for n, fd in fds if fd.write]
    if writable:
        L += ["   template< typename WriterType >",
              "   void writeVariables( WriterType& writer, const GlobalIndexType& numberOfParticles, const GlobalIndexType firstActiveParticle = 0 )",
              "   {"]
        for n, fd in writable:
            cap = _vtk_name(n)
            if fd.dtype == "scalar":
                L.append(f'      writer.template writePointData< ScalarArrayType >( {n}, "{cap}", numberOfParticles, firstActiveParticle, 1 );')
            else:
                L.append(f'      writer.template writeVector< VectorArrayType, RealType >( {n}, "{cap}", numberOfParticles, firstActiveParticle, 3 );')
        L += ["   }", ""]
    # MPI synchronization
    syncable = [(n, fd) for n, fd in fds if fd.t >= 2]
    if syncable:
        L += ["#ifdef HAVE_MPI",
              "   template< typename Synchronizer, typename DistributedParticlesPointer >",
              "   void synchronizeVariables( Synchronizer& synchronizer,",
              "                             DistributedParticlesPointer& distributedParticles )",
              "   {"]
        for n, _ in syncable:
            L.append(f"      synchronizer.synchronize( {n}, distributedParticles );")
        L += ["   }", "#endif", ""]
    L.append("};")
    return L

def _boundary_var_class(boundary_cls):
    """Generate BoundaryVariables that inherits from FluidVariables."""
    # For DBC, boundary inherits all fluid variables (no specialization needed)
    return ["template< typename SPHState >",
            "class BoundaryVariables : public FluidVariables< SPHState >",
            "{};"]

def _vtk_name(field_name):
    """Map field name to VTK variable name (matching hand-written convention)."""
    names = {"rho": "Density", "v": "Velocity", "p": "Pressure",
             "drho": "DensityRate", "a": "Acceleration"}
    return names.get(field_name, field_name.capitalize())

# ─── control.h ───────────────────────────────────────────────────────────────
def _gen_control(constants_cls, namespace):
    sc = constants_cls.scalar_names()
    vc = constants_cls.vector_names()
    # Add derived constants needed by the framework
    extra_scalars = [c for c in ["coefB", "p0"] if c not in sc]
    all_sc = sc + extra_scalars

    L = ["// control.h  –  AUTO-GENERATED", "#pragma once", "",
         "#include <SPH/SPHTraits.h>",
         "#include <SPH/TimeStep.h>",
         "#include <SPH/Kernels.h>",
         "#include <SPH/Models/EquationOfState.h>",
         "#include <SPH/Models/DiffusiveTerms.h>",
         "#include <SPH/Models/VisousTerms.h>",
         "#include <SPH/Models/WCSPH_DBC/BoundaryConditionsTypes.h>",
         "#include <SPH/Models/WCSPH_DBC/IntegrationSchemes/VerletScheme.h>",
         '#include <SPH/Models/WCSPH_DBC/OpenBoundaryConfig.h>',
         "",
         "namespace TNL {", "namespace SPH {", "",
         f"// Parameters for {namespace}",
         "template< typename SPHDefs >",
         "class WCSPH_DBCConfig", "{", "public:",
         "   using SPHConfig         = typename SPHDefs::SPHConfig;",
         "   using SPHTraitsType     = SPHFluidTraits< SPHConfig >;",
         "   using RealType          = typename SPHTraitsType::RealType;",
         "   using VectorType        = typename SPHTraitsType::VectorType;",
         "   using KernelFunction    = typename SPHDefs::KernelFunction;",
         "   using DiffusiveTerm     = typename SPHDefs::DiffusiveTerm;",
         "   using ViscousTerm       = typename SPHDefs::ViscousTerm;",
         "   using EOS               = typename SPHDefs::EOS;",
         "   using BCType            = typename SPHDefs::BCType;",
         "   using IntegrationScheme = typename SPHDefs::IntegrationScheme;",
         "   using TimeStepping      = typename SPHDefs::TimeStepping;",
         ""]

    # configSetupModel
    L += ["   static void configSetupModel( TNL::Config::ConfigDescription& config )",
          "   {",
          f'      config.addDelimiter( "{namespace} model parameters" );']
    for n in sc:
        desc = _param_desc(n)
        L.append(f'      config.addEntry< double >( "{n}", "{desc}", 0 );')
    for n in vc:
        for ax in ("x", "y", "z"):
            L.append(f'      config.addEntry< double >( "{n}-{ax}", "{_param_desc(n)} {ax}", 0 );')
    L += ["   }", ""]

    # init
    L += ["   void init( TNL::Config::ParameterContainer& parameters )", "   {"]
    for n in sc:
        L.append(f'      {n} = parameters.getParameter< double >( "{n}" );')
    for n in vc:
        L.append(f'      {n} = parameters.getXyz< VectorType >( "{n}" );')
    # Derived constants
    if "coefB" not in sc:
        L.append("      coefB = speedOfSound * speedOfSound * rho0 / 7.f;")
    if "dtMin" in sc and "dtInit" in sc:
        L.append("      dtMin = 0.05f * dtInit;")
    L += ["   }", ""]

    # Member declarations
    for n in all_sc:
        L.append(f"   RealType   {n} = 0.f;")
    for n in vc:
        L.append(f"   VectorType {n} = 0.f;")
    L += ["", "};", ""]

    # writePrologModel — free function called by SPHMultiset_CFD via ADL
    L += _gen_write_prolog(sc, vc, namespace)

    L += ["",
          "} // SPH", "} // TNL", "", "// end control.h"]
    return "\n".join(L)


def _gen_write_prolog(sc, vc, namespace):
    L = [
        "template< typename ModelParams >",
        "void",
        "writePrologModel( TNL::Logger& logger, ModelParams& modelParams )",
        "{",
        f'   logger.writeHeader( "TNL::SPH::{namespace} model parameters" );',
        '   logger.writeParameter( "Resolution parameters", "" );',
    ]
    prolog_scalars = ["dp", "h", "mass"]
    for n in prolog_scalars:
        if n in sc:
            desc = _prolog_label(n)
            L.append(f'   logger.writeParameter( "{desc}:", modelParams.{n}, 1 );')
    if "h" in sc and "dp" in sc:
        L.append('   logger.writeParameter( "Spatial resolution (h/dp):", modelParams.h / modelParams.dp, 1 );')

    L.append('   logger.writeParameter( "Model parameters", "" );')
    remaining_sc = [n for n in sc if n not in prolog_scalars]
    for n in remaining_sc:
        desc = _prolog_label(n)
        L.append(f'   logger.writeParameter( "{desc}:", modelParams.{n}, 1 );')
    for n in vc:
        desc = _prolog_label(n)
        L.append(f'   logger.writeParameter( "{desc}:", modelParams.{n} );')

    L += ['   logger.writeParameter( "External bulk force:", modelParams.gravity );',
          "}"]
    return L


def _prolog_label(name):
    labels = {
        "h": "Smoothing length (h)",
        "dp": "Initial particle distance (dp)",
        "mass": "Particle mass (mass)",
        "delta": "Diffusive term coefficient (delta)",
        "alpha": "Artificial viscosity coefficient (alpha)",
        "dynamicViscosity": "Dynamic viscosity (dynamicViscosity)",
        "speedOfSound": "Speed of sound (speedOfSound)",
        "rho0": "Referential density (rho0)",
        "p0": "Background pressure (p0)",
        "dtInit": "Initial time step (dtInit)",
        "cfl": "CFL number (CFL)",
        "dtMin": "Minimal time step (dtMin)",
        "eps": "Epsilon (eps)",
        "gravity": "External bulk force",
        "coefB": "Coefficient of EOS (coefB)",
    }
    return labels.get(name, name)

def _param_desc(name):
    descs = {
        "h": "SPH method smoothing length.",
        "dp": "Initial particle distance.",
        "mass": "Mass of particle, constant for all particles.",
        "delta": "Coefficient of artificial delta-WCSPH diffusive term.",
        "alpha": "Coefficient of artificial viscous term.",
        "dynamicViscosity": "Dynamic viscosity coefficient.",
        "speedOfSound": "Numerical speed of sound.",
        "rho0": "Referential density of the medium.",
        "p0": "Background pressure.",
        "dtInit": "Initial time step.",
        "cfl": "CFL number.",
        "dtMin": "Minimal allowed time step.",
        "eps": "Coefficient to prevent denominator from zero.",
        "gravity": "External bulk forces.",
    }
    return descs.get(name, "")

# ─── Interactions.h ─────────────────────────────────────────────────────────
def _gen_interactions_h(interactions_cls, namespace):
    L = ["// Interactions.h  –  AUTO-GENERATED", "#pragma once", "",
         '#include "../../SPHTraits.h"',
         '#include "../WCSPH_DBC/BoundaryConditionsTypes.h"',
         '#include "../WCSPH_DBC/OpenBoundaryConfig.h"',
         '#include "../WCSPH_DBC/OpenBoundaryConditions.h"',
         '#include "Variables.h"',
         '#include <TNL/Matrices/StaticMatrix.h>',
         '#include "../EquationOfState.h"',
         '#include "../DiffusiveTerms.h"',
         '#include "../VisousTerms.h"',
         '#include "control.h"',
         "",
         "namespace TNL {", "namespace SPH {", "",
         "template< typename Particles, typename ModelConfig >",
         f"class {namespace}", "{", "public:",
         f"   using Model = {namespace}< Particles, ModelConfig >;",
         "   using ModelParams = WCSPH_DBCConfig< ModelConfig >;",
         "   using ParticlesType = Particles;",
         "   using ModelConfigType = ModelConfig;",
         "   using SPHConfig = typename ModelConfig::SPHConfig;",
         "   using DeviceType = typename SPHConfig::DeviceType;",
         "",
         "   using SPHTraitsType = SPHFluidTraits< SPHConfig >;",
         "   using LocalIndexType = typename SPHTraitsType::LocalIndexType;",
         "   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;",
         "   using IndexVectorType = typename SPHTraitsType::IndexVectorType;",
         "   using RealType = typename SPHTraitsType::RealType;",
         "   using VectorType = typename SPHTraitsType::VectorType;",
         "   using Matrix = Matrices::StaticMatrix< RealType, SPHConfig::spaceDimension + 1, SPHConfig::spaceDimension + 1 >;",
         "   using VectorExtendedType = Containers::StaticVector< SPHConfig::spaceDimension + 1, RealType >;",
         "",
         "   using FluidVariables = FluidVariables< SPHConfig >;",
         "   using BoundaryVariables = BoundaryVariables< ModelConfig >;",
         "   using OpenBoundaryVariables = OpenBoundaryVariables< SPHConfig >;",
         "   using IntegrationSchemeType = typename ModelConfig::IntegrationScheme;",
         "   using IntegrationSchemeVariables = typename IntegrationSchemeType::IntegrationSchemeVariablesType;",
         "   using KernelFunction = typename ModelConfig::KernelFunction;",
         "   using DiffusiveTerm = typename ModelConfig::DiffusiveTerm;",
         "   using ViscousTerm = typename ModelConfig::ViscousTerm;",
         "   using EOS = typename ModelConfig::EOS;",
         "",
         "   using OpenBoundaryConfig = DBCOpenBoundaryConfig< SPHConfig >;",
         "   using OpenBoundaryModel = OpenBoundaryConditionsBuffers< SPHConfig, ModelConfig >;",
         "",
         f'   static std::string writeModelType() {{ return "TNL::SPH::{namespace}"; }}',
         f"   {namespace}() = default;", "",
         # Method declarations
         "   template< typename EOS_ = EOS, typename PhysicalObjectPointer >",
         "   void computePressureFromDensity( PhysicalObjectPointer& physicalObject, ModelParams& modelParams );",
         "",
         "   template< typename FluidPointer, typename BoundaryPointer >",
         "   void interaction( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams );",
         "",
         "   template< typename FluidPointer, typename BoundaryPointer >",
         "   void updateSolidBoundary( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams );",
         "",
         "   template< typename FluidPointer, typename BoundaryPointer >",
         "   void initializeInteraction( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams ) {}",
         "",
         "   template< typename FluidPointer, typename BoundaryPointer >",
         "   void finalizeInteraction( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams ) {}",
         "",
         "   template< typename FluidPointer, typename BoundaryPointer >",
         "   void initializeBoundaryInteraction( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams ) {}",
         "",
         "   template< typename FluidPointer, typename BoundaryPointer >",
         "   void finalizeBoundaryInteraction( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams ) {}",
         "",
         # Open boundary stubs (needed for compilation, never called without open boundary patches)
         "   template< typename OpenBoundaryPointer, typename BoundaryPointer >",
         "   void updateSolidBoundaryOpenBoundary( BoundaryPointer& boundary, OpenBoundaryPointer& openBoundary, ModelParams& modelParams ) {}",
         "",
         "   template< typename FluidPointer, typename OpenBoundaryPointer >",
         "   void interactionWithOpenBoundary( FluidPointer& fluid, OpenBoundaryPointer& openBoundary, ModelParams& modelParams ) {}",
         "",
         "   template< typename FluidPointer, typename OpenBoundaryPointer >",
         "   void interactionWithBoundaryPatches( FluidPointer& fluid, OpenBoundaryPointer& boundaryPatch, ModelParams& modelParams ) {}",
         "",
         "   template< typename FluidPointer, typename OpenBoundaryPointer >",
         "   void extrapolateOpenBoundaryData( FluidPointer& fluid, OpenBoundaryPointer& openBoundary, ModelParams& modelParams, OpenBoundaryConfig& openBoundaryParams ) {}",
         "",
         "};", "",
         "} // SPH", "} // TNL", "",
         '#include "Interactions.hpp"',
         "", "// end Interactions.h"]
    return "\n".join(L)

# ─── Interactions.hpp ────────────────────────────────────────────────────────
def _gen_interactions_hpp(interactions_cls, fluid_cls, boundary_cls, constants_cls, namespace):
    L = ["// Interactions.hpp  –  AUTO-GENERATED",
         '#include "Interactions.h"', "",
         "namespace TNL {", "namespace SPH {", ""]

    phases = _phase_map(interactions_cls)
    for pname, pairs in phases.items():
        L += _gen_phase_method(pname, pairs, fluid_cls, boundary_cls, constants_cls, namespace)
        L.append("")

    # computePressureFromDensity
    L += _gen_compute_pressure(namespace)
    L.append("")

    L += ["} // SPH", "} // TNL", "", "// end Interactions.hpp"]
    return "\n".join(L)

def _gen_compute_pressure(namespace):
    return [
        "template< typename Particles, typename ModelConfig >",
        "template< typename EOS_, typename PhysicalObjectPointer >",
        "void",
        f"{namespace}< Particles, ModelConfig >::computePressureFromDensity(",
        "   PhysicalObjectPointer& physicalObject, ModelParams& modelParams )",
        "{",
        "   auto view_rho = physicalObject->getVariables()->rho.getView();",
        "   auto view_p = physicalObject->getVariables()->p.getView();",
        "   typename EOS_::ParamsType eosParams( modelParams );",
        "",
        "   auto evalPressure = [=] __cuda_callable__ ( int i ) mutable",
        "   {",
        "      view_p[ i ] = EOS_::DensityToPressure( view_rho[ i ], eosParams );",
        "   };",
        "   physicalObject->getParticles()->forAll( evalPressure );",
        "}"]

def _gen_phase_method(pname, pairs, fluid_cls, boundary_cls, constants_cls, namespace):
    """Generate one phase method (interaction or updateSolidBoundary)."""
    owner = pairs[0].set_i  # "fluid" or "boundary"
    o_obj = "fluid" if owner == "fluid" else "boundary"
    o_sfx = _view_suffix(owner)

    L = ["template< typename Particles, typename ModelConfig >",
         "template< typename FluidPointer, typename BoundaryPointer >",
         "void",
         f"{namespace}< Particles, ModelConfig >::{pname}(",
         "   FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams )",
         "{"]

    # ─── Search tokens ──────────────────────────────────────────────────────
    seen_tokens = set()
    for p in pairs:
        tok = _search_token(owner, p.set_j)
        if tok not in seen_tokens:
            seen_tokens.add(tok)
            tgt_obj = "fluid" if p.set_j == "fluid" else "boundary"
            L.append(f"   auto {tok} = {o_obj}->getParticles()->getSearchToken( {tgt_obj}->getParticles() );")
    L.append("")

    # ─── Constants ──────────────────────────────────────────────────────────
    L.append("   const RealType searchRadius = fluid->getParticles()->getSearchRadius();")
    if constants_cls:
        for n in constants_cls.scalar_names():
            L.append(f"   const RealType {n} = modelParams.{n};")
        for n in constants_cls.vector_names():
            L.append(f"   const VectorType {n} = modelParams.{n};")
    L.append("   typename EOS::ParamsType eosParams( modelParams );")
    L.append("")

    # ─── Views ──────────────────────────────────────────────────────────────
    written_fields = set()
    for p in pairs:
        for s in p.stmts:
            written_fields.add(s.field)

    j_accessed_fields = set()
    for p in pairs:
        j_accessed_fields |= _j_fields(p, fluid_cls, boundary_cls)

    all_sets = set()
    for p in pairs:
        all_sets.add(p.set_i)
        all_sets.add(p.set_j)

    seen_views = set()
    def add_view(field_name, pset, mutable):
        key = (field_name, pset)
        if key in seen_views:
            return
        seen_views.add(key)
        sfx = _view_suffix(pset)
        obj = "fluid" if pset == "fluid" else "boundary"
        q = "auto" if mutable else "const auto"
        m = "getView()" if mutable else "getConstView()"
        if field_name == "points":
            L.append(f"   {q} view_points{sfx} = {obj}->getParticles()->getPoints().{m};")
        else:
            display = field_name[0].upper() + field_name[1:] if mutable else field_name
            L.append(f"   {q} view_{display}{sfx} = {obj}->getVariables()->{field_name}.{m};")

    for pset in all_sets:
        cls = fluid_cls if pset == "fluid" else boundary_cls

        add_view("points", pset, False)

        if pset == owner:
            for n, fd in _fields(cls):
                is_written = n in written_fields
                if fd.eos and not is_written:
                    continue
                if is_written:
                    add_view(n, pset, True)
                elif n in ("v", "rho"):
                    add_view(n, pset, False)
        else:
            for n, fd in _fields(cls):
                if fd.eos:
                    continue
                if n in j_accessed_fields:
                    add_view(n, pset, False)

    L.append("")

    # ─── Lambdas ────────────────────────────────────────────────────────────
    for p in pairs:
        L += _gen_lambda(p, fluid_cls, boundary_cls)
        L.append("")

    # ─── Particle loop ──────────────────────────────────────────────────────
    L += _gen_particle_loop(pairs, owner, o_sfx, o_obj, fluid_cls, boundary_cls, constants_cls)

    L.append("}")
    return L

def _gen_lambda(pair, fluid_cls, boundary_cls):
    """Generate one neighbor-loop lambda."""
    j_sfx = _view_suffix(pair.set_j)
    outs = _output_fields(pair)

    # Determine field types
    all_fields = {}
    for cls in [fluid_cls, boundary_cls]:
        for n, fd in _fields(cls):
            all_fields[n] = fd

    # Lambda parameters: r_i, v_i, rho_i, p_i (inputs), then output pointers
    input_params = ["VectorType& r_i", "VectorType& v_i", "RealType& rho_i", "RealType& p_i"]
    output_params = []
    for f in outs:
        fd = all_fields.get(f)
        ct = _ctype_short(fd) if fd else "RealType"
        output_params.append(f"{ct}* {f}_i")

    params_str = ",\n         ".join(
        ["LocalIndexType i, LocalIndexType j"] + input_params + output_params
    )

    L = [f"   auto {pair.name} = [=] __cuda_callable__ (",
         f"         {params_str} ) mutable",
         "   {",
         f"      const VectorType r_j = view_points{j_sfx}[ j ];",
         "      const VectorType r_ij = r_i - r_j;",
         "      const RealType drs = l2Norm( r_ij );",
         "      if( drs <= searchRadius )",
         "      {"]

    # Load j-particle data from views
    j_cls = fluid_cls if pair.set_j == "fluid" else boundary_cls
    j_fields_used = _j_fields(pair, fluid_cls, boundary_cls)

    for n, fd in _fields(j_cls):
        if n in j_fields_used or n in ("v", "rho"):  # Always load v_j and rho_j
            if fd.eos:
                continue  # EOS fields computed below
            ct = _ctype_short(fd)
            L.append(f"         const {ct} {n}_j = view_{n}{j_sfx}[ j ];")

    # Compute p_j from EOS (always, matching hand-written code)
    L.append("         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );")

    # Kernel gradient
    L += ["         const RealType F = KernelFunction::F( drs, h );",
          "         const VectorType gradW = r_ij * F;", ""]

    # Locals (let declarations)
    for d in pair.locals:
        L.append(f"         const auto {d.name} = {emit(d.expr)};")
    if pair.locals:
        L.append("")

    # Accumulation statements
    for s in pair.stmts:
        op = "+=" if s.kind == "add" else "="
        L.append(f"         *{s.field}_i {op} {emit(s.expr)};")

    L += ["      }", "   };"]
    return L

def _gen_particle_loop(pairs, owner, o_sfx, o_obj, fluid_cls, boundary_cls, constants_cls):
    """Generate the main particle loop."""
    # Collect all output fields across all pairs
    all_outs = []
    for p in pairs:
        for f in _output_fields(p):
            if f not in all_outs:
                all_outs.append(f)

    # Field types
    all_fields = {}
    for cls in [fluid_cls, boundary_cls]:
        for n, fd in _fields(cls):
            all_fields[n] = fd

    has_gravity = constants_cls and "gravity" in constants_cls.vector_names()

    L = ["   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable",
         "   {",
         f"      const VectorType r_i = view_points{o_sfx}[ i ];",
         f"      const VectorType v_i = view_v{o_sfx}[ i ];",
         f"      const RealType rho_i = view_rho{o_sfx}[ i ];",
         "      const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );",
         ""]

    # Initialize output accumulators
    for f in all_outs:
        fd = all_fields.get(f)
        z = _init_val(fd) if fd else "0.f"
        L.append(f"      {_ctype_short(fd) if fd else 'RealType'} {f}_i = {z};")
    L.append("")

    # Neighbor loops
    for p in pairs:
        ptrs = ", ".join(f"&{f}_i" for f in _output_fields(p))
        tok = _search_token(owner, p.set_j)
        loop = "Particles::NeighborsLoop" if p.set_j == owner else "Particles::NeighborsLoopAnotherSet"
        L.append(f"      {loop}::exec( i, r_i, {tok}, {p.name},")
        L.append(f"                        v_i, rho_i, p_i, {ptrs} );")

    L.append("")
    # Write outputs (with gravity addition for acceleration field)
    for f in all_outs:
        fd = all_fields.get(f)
        if f == "a" and has_gravity:
            L.append(f"      {f}_i += gravity;")
        # Use capitalized name for mutable views
        display = f[0].upper() + f[1:] if f in all_outs else f
        L.append(f"      view_{display}{o_sfx}[ i ] = {f}_i;")

    L += ["   };",
          f"   {o_obj}->getParticles()->forAll( particleLoop );"]

    # Periodic patches
    L += ["",
          f"   if( {o_obj}->periodicPatches.size() > 0 ) {{",
          "      for( long unsigned int p = 0; p < std::size( " + o_obj + "->periodicPatches ); p++ ) {",
          "",
          f"         const auto zoneParticleIndices_view = {o_obj}->periodicPatches[ p ]->particleZone.getParticlesInZone().getConstView();",
          f"         const GlobalIndexType numberOfZoneParticles = {o_obj}->periodicPatches[ p ]->particleZone.getNumberOfParticles();",
          f"         const VectorType shift = {o_obj}->periodicPatches[ p ]->config.shift;",
          "",
          "         auto periodicParticleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable",
          "         {",
          "            const GlobalIndexType pi = zoneParticleIndices_view[ i ];",
          f"            const VectorType r_i = view_points{o_sfx}[ pi ] + shift;",
          f"            const VectorType v_i = view_v{o_sfx}[ pi ];",
          f"            const RealType rho_i = view_rho{o_sfx}[ pi ];",
          "            const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );",
          ""]

    # Initialize accumulators for periodic
    for f in all_outs:
        fd = all_fields.get(f)
        z = _init_val(fd) if fd else "0.f"
        L.append(f"            {_ctype_short(fd) if fd else 'RealType'} {f}_i = {z};")
    L.append("")

    # Neighbor loops for periodic
    for p in pairs:
        ptrs = ", ".join(f"&{f}_i" for f in _output_fields(p))
        tok = _search_token(owner, p.set_j)
        loop = "Particles::NeighborsLoop" if p.set_j == owner else "Particles::NeighborsLoopAnotherSet"
        L.append(f"            {loop}::exec( pi, r_i, {tok}, {p.name},")
        L.append(f"                              v_i, rho_i, p_i, {ptrs} );")

    L.append("")
    for f in all_outs:
        display = f[0].upper() + f[1:]
        L.append(f"            view_{display}{o_sfx}[ pi ] += {f}_i;")

    L += ["         };",
          "         Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, periodicParticleLoop );",
          "      }",
          "   }"]

    return L
