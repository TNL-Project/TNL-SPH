"""
sph_codegen/scheme.py
=====================
Field descriptors, particle-set metaclasses, @interactions / @integration decorators.

Both decorators parse function bodies via Python AST — no execution, no
monkey-patching, no namespace injection. Pure structural analysis.
"""
from __future__ import annotations
import ast
import inspect
import textwrap
from .dsl import (
    Expr, Lit, _B, _U, _C, Fn, Dot, Ternary,
    ConstRef, FieldRef, LocalRef, TimeRef,
    Stmt, LocalDecl, _w,
)

# -- Field descriptors -------------------------------------------------------

class _FieldDesc:
    def __init__(self, dtype, t=1, read=False, write=False, eos=False):
        self.dtype  = dtype
        self.t      = t
        self.read   = read
        self.write  = write
        self.eos    = eos
        self.name   = None
        self.pset   = None

def ScalarField(t=1, read=False, write=False, eos=False): return _FieldDesc("scalar", t, read, write, eos)
def VectorField(t=1, read=False, write=False, eos=False): return _FieldDesc("vector", t, read, write, eos)

# -- ParticleSet metaclass ---------------------------------------------------

class _PSMeta(type):
    def __new__(mcs, name, bases, ns):
        fields = {k: v for k, v in ns.items() if isinstance(v, _FieldDesc)}
        cls = super().__new__(mcs, name, bases, ns)
        cls._fields = fields
        return cls

class FluidSet(metaclass=_PSMeta):
    _pset = "fluid"

class BoundarySet(metaclass=_PSMeta):
    _pset = "boundary"

# -- Constants ---------------------------------------------------------------

class ConstantSet:
    Scalars = ""
    Vectors = ""
    @classmethod
    def scalar_names(cls): return cls.Scalars.split() if cls.Scalars else []
    @classmethod
    def vector_names(cls): return cls.Vectors.split() if cls.Vectors else []
    @classmethod
    def all_names(cls): return cls.scalar_names() + cls.vector_names()

# -- Index singletons (for signature inspection) -----------------------------

class _Idx:
    def __init__(self, n): self._n = n
    def __repr__(self): return self._n

i = _Idx("i")
j = _Idx("j")

# -- Helpers exported for LSP only (never executed by the parser) ------------

def dot(a, b):  return Dot(_w(a), _w(b))
def norm(v):    return Fn("l2Norm", [_w(v)])
def powf(b, e): return Fn("powf", [_w(b), _w(e)])

gradW = LocalRef("gradW")
r_ij  = LocalRef("r_ij")
drs   = LocalRef("drs")
t     = LocalRef("t")
dt    = LocalRef("dt")

# -- AST expression builder --------------------------------------------------

_BINOP_MAP = {
    ast.Add: "+", ast.Sub: "-", ast.Mult: "*", ast.Div: "/",
}
_CMPOP_MAP = {
    ast.Lt: "<", ast.LtE: "<=", ast.Gt: ">", ast.GtE: ">=",
    ast.Eq: "==", ast.NotEq: "!=",
}

_BUILTINS = {"gradW", "r_ij", "drs", "searchRadius", "t", "dt", "step"}


class _ExprBuilder:
    """Walks Python AST expression nodes and produces Expr trees."""

    def __init__(self, field_names, const_names, si, sj):
        self.field_names = field_names
        self.const_names = const_names
        self.si = si
        self.sj = sj
        self.known_locals: set[str] = set()

    def build(self, node) -> Expr:
        if isinstance(node, ast.Constant):
            v = node.value
            return Lit(float(v) if isinstance(v, int) else v)

        if isinstance(node, ast.Name):
            return self._build_name(node.id)

        if isinstance(node, ast.Subscript):
            return self._build_subscript(node)

        if isinstance(node, ast.BinOp):
            op = _BINOP_MAP.get(type(node.op))
            if op is None:
                raise ValueError(f"Unsupported binary operator: {type(node.op).__name__}")
            return _B(op, self.build(node.left), self.build(node.right))

        if isinstance(node, ast.UnaryOp):
            if isinstance(node.op, ast.USub):
                return _U("-", self.build(node.operand))
            if isinstance(node.op, ast.UAdd):
                return self.build(node.operand)
            raise ValueError(f"Unsupported unary operator: {type(node.op).__name__}")

        if isinstance(node, ast.Compare):
            if len(node.ops) != 1:
                raise ValueError("Chained comparisons not supported")
            op = _CMPOP_MAP.get(type(node.ops[0]))
            if op is None:
                raise ValueError(f"Unsupported comparison: {type(node.ops[0]).__name__}")
            return _C(op, self.build(node.left), self.build(node.comparators[0]))

        if isinstance(node, ast.IfExp):
            return Ternary(
                self.build(node.test),
                self.build(node.body),
                self.build(node.orelse),
            )

        if isinstance(node, ast.Call):
            return self._build_call(node)

        raise ValueError(f"Unsupported AST node: {type(node).__name__}")

    def _build_name(self, name_id: str) -> Expr:
        if name_id in self.known_locals:
            return LocalRef(name_id)
        if name_id in _BUILTINS:
            return LocalRef(name_id)
        if name_id in self.const_names:
            return ConstRef(name_id)
        raise ValueError(f"Unknown name '{name_id}' in expression")

    def _build_subscript(self, node) -> Expr:
        if not isinstance(node.value, ast.Name):
            raise ValueError("Subscript target must be a field name")
        field_name = node.value.id
        if field_name not in self.field_names:
            raise ValueError(f"'{field_name}' is not a declared field")
        idx = _get_slice_id(node.slice)
        if idx not in ("i", "j"):
            raise ValueError(f"Field index must be 'i' or 'j', got '{idx}'")
        pset = self.si if idx == "i" else self.sj
        return FieldRef(field_name, idx, pset)

    def _build_call(self, node) -> Expr:
        if isinstance(node.func, ast.Subscript):
            return self._build_time_ref(node)

        if not isinstance(node.func, ast.Name):
            raise ValueError("Only simple function calls supported")
        fname = node.func.id
        args = [self.build(a) for a in node.args]

        if fname == "dot":
            if len(args) != 2:
                raise ValueError("dot() requires exactly 2 arguments")
            return Dot(args[0], args[1])
        if fname == "norm":
            if len(args) != 1:
                raise ValueError("norm() requires exactly 1 argument")
            return Fn("l2Norm", args)
        if fname == "powf":
            if len(args) != 2:
                raise ValueError("powf() requires exactly 2 arguments")
            return Fn("powf", args)
        raise ValueError(f"Unknown function '{fname}'")

    def _build_time_ref(self, node) -> Expr:
        """field[idx](t+dt) — integration time reference on RHS."""
        sub = node.func
        field_name = sub.value.id
        idx = _get_slice_id(sub.slice)
        pset = self.si if idx == "i" else self.sj
        offset = _parse_t_offset(node.args[0])
        return TimeRef(field_name, pset, offset)


def _get_slice_id(slice_node) -> str:
    if isinstance(slice_node, ast.Name):
        return slice_node.id
    if hasattr(ast, "Index") and isinstance(slice_node, ast.Index):
        return slice_node.value.id
    raise ValueError("Field index must be a simple name (i or j)")


def _parse_t_offset(node) -> int:
    if isinstance(node, ast.Name):
        return 0
    if isinstance(node, ast.BinOp):
        if isinstance(node.op, ast.Add):
            return +1
        if isinstance(node.op, ast.Sub):
            return -1
    return 0


# -- Interaction records -----------------------------------------------------

class _IRecord:
    def __init__(self, name, si, sj, locals_decls, stmts):
        self.name = name
        self.set_i = si
        self.set_j = sj
        self.locals = locals_decls
        self.stmts = stmts


_PSET_MAP = {"fluid": "fluid", "boundary": "boundary", "bound": "boundary"}


def _detect_sets(fname):
    parts = fname.split("_")
    found = [_PSET_MAP[p] for p in parts if p in _PSET_MAP]
    si = found[0] if len(found) > 0 else "fluid"
    sj = found[1] if len(found) > 1 else "fluid"
    return si, sj


def _field_names(cls):
    return [k for k, v in vars(cls).items() if isinstance(v, _FieldDesc)]


def _parse_interaction(fn, fluid_cls, boundary_cls, constants_cls, si, sj):
    source = textwrap.dedent(inspect.getsource(fn))
    tree = ast.parse(source)
    func_def = tree.body[0]

    field_names = set(_field_names(fluid_cls)) | set(_field_names(boundary_cls))
    const_names = set(constants_cls.all_names()) if constants_cls else set()

    builder = _ExprBuilder(field_names, const_names, si, sj)
    locals_decls: list[LocalDecl] = []
    stmts: list[Stmt] = []

    for node in func_def.body:
        _process_interaction_stmt(node, builder, locals_decls, stmts)

    return _IRecord(fn.__name__, si, sj, locals_decls, stmts)


def _process_interaction_stmt(node, builder, locals_decls, stmts):
    if isinstance(node, ast.Assign):
        if len(node.targets) != 1:
            raise ValueError("Multiple assignment targets not supported")
        target = node.targets[0]
        if not isinstance(target, ast.Name):
            raise ValueError("Assignment target must be a simple name")
        name = target.id
        expr = builder.build(node.value)
        builder.known_locals.add(name)
        locals_decls.append(LocalDecl(name, expr))

    elif isinstance(node, ast.AugAssign):
        if not isinstance(node.op, ast.Add):
            raise ValueError(f"Only += supported, got {type(node.op).__name__}=")
        target = node.target
        if not isinstance(target, ast.Subscript):
            raise ValueError("+= target must be field[index]")
        field_name = target.value.id
        idx = _get_slice_id(target.slice)
        pset = builder.si if idx == "i" else builder.sj
        expr = builder.build(node.value)
        stmts.append(Stmt("add", field_name, idx, pset, expr))

    elif isinstance(node, ast.If):
        for s in node.body:
            _process_interaction_stmt(s, builder, locals_decls, stmts)
        for s in node.orelse:
            _process_interaction_stmt(s, builder, locals_decls, stmts)

    elif isinstance(node, ast.Pass):
        pass

    else:
        raise ValueError(f"Unsupported statement: {type(node).__name__}")


def interactions(fluid_cls, boundary_cls, constants_cls=None):
    def decorator(cls):
        records = []
        for fname, fn in vars(cls).items():
            if fname.startswith("_") or not callable(fn):
                continue
            sig = inspect.signature(fn)
            params = list(sig.parameters.keys())
            if set(params) - {"i", "j"}:
                continue
            si, sj = _detect_sets(fname)
            records.append(_parse_interaction(fn, fluid_cls, boundary_cls, constants_cls, si, sj))
        cls._interactions = records
        return cls
    return decorator


# -- Integration records -----------------------------------------------------

class _IntStmt:
    def __init__(self, field, pset, offset, rhs_expr, condition=None):
        self.field = field
        self.pset = pset
        self.offset = offset
        self.rhs_expr = rhs_expr
        self.condition = condition


class _IntRecord:
    def __init__(self, fname, pset, stmts):
        self.name = fname
        self.pset = pset
        self.stmts = stmts


def _parse_integration(fn, fluid_cls, boundary_cls, constants_cls, pset):
    source = textwrap.dedent(inspect.getsource(fn))
    tree = ast.parse(source)
    func_def = tree.body[0]

    field_names = set(_field_names(fluid_cls)) | set(_field_names(boundary_cls))
    const_names = set(constants_cls.all_names()) if constants_cls else set()

    builder = _ExprBuilder(field_names, const_names, pset, pset)
    stmts: list[_IntStmt] = []

    def process_body(nodes, cond_str):
        for node in nodes:
            if isinstance(node, ast.If):
                cond = ast.unparse(node.test)
                process_body(node.body, cond)
                process_body(node.orelse, None)
            elif isinstance(node, ast.Expr):
                val = node.value
                if isinstance(val, ast.BinOp) and isinstance(val.op, ast.LShift):
                    _parse_integration_lshift(val, builder, pset, stmts, cond_str)
            elif isinstance(node, ast.Pass):
                pass

    process_body(func_def.body, None)
    return _IntRecord(fn.__name__, pset, stmts)


def _parse_integration_lshift(node, builder, pset, stmts, cond_str):
    tgt = node.left
    rhs = node.right

    if not isinstance(tgt, ast.Call) or not isinstance(tgt.func, ast.Subscript):
        raise ValueError("Integration target must be field[idx](t+dt)")

    field_name = tgt.func.value.id
    offset = _parse_t_offset(tgt.args[0])
    rhs_expr = builder.build(rhs)
    stmts.append(_IntStmt(field_name, pset, offset, rhs_expr, cond_str))


def integration(fluid_cls, boundary_cls, constants_cls=None):
    def decorator(cls):
        records = []
        for fname, fn in vars(cls).items():
            if fname.startswith("_") or not callable(fn):
                continue
            sig = inspect.signature(fn)
            params = list(sig.parameters.keys())
            if params not in [["i"], ["self", "i"]]:
                continue
            pset = "boundary" if "boundary" in fname else "fluid"
            records.append(_parse_integration(fn, fluid_cls, boundary_cls, constants_cls, pset))
        cls._rules = records
        return cls
    return decorator
