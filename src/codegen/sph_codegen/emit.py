"""
sph_codegen/emit.py
====================
Walks dsl Expr trees → C++ expression strings.

Key principle: FieldRef nodes emit as local lambda-parameter names
(e.g.  rho_i, v_j) — NOT as raw view accesses. The generator is
responsible for declaring these locals and loading them from views.
"""
from __future__ import annotations
from .dsl import (
    Expr, Lit, _B, _U, _C, Fn, Dot, Ternary,
    ConstRef, FieldRef, LocalRef, TimeRef,
    Stmt, LocalDecl,
)

_PREC = {
    "||":1,"&&":2,"==":3,"!=":3,"<":3,"<=":3,">":3,">=":3,
    "+":4,"-":4,"*":5,"/":5,"u-":7,
}

def _p(node):
    if isinstance(node, _B): return _PREC.get(node.op,5)
    if isinstance(node, _C): return _PREC.get(node.op,3)
    if isinstance(node, _U): return _PREC.get("u-",7)
    return 99

def emit(node, pp=0) -> str:
    if node is None: return "0.f"

    # ── Literal ──────────────────────────────────────────────────────────────
    if isinstance(node, Lit):
        v = node.value
        if isinstance(v, float):
            if v == int(v):
                return f"{int(v)}.f"
            return f"{v}f"
        return str(v)

    # ── Constant reference (modelParams member) ──────────────────────────────
    if isinstance(node, ConstRef):
        return node.name

    # ── Local variable (let-declared inside interaction) ─────────────────────
    if isinstance(node, LocalRef):
        return node.name

    # ── Field reference: rho(i) → rho_i, v(j) → v_j ──────────────────────────
    #  These are lambda parameters (for i) or loaded locals (for j).
    #  The generator declares them; the emitter just uses the name.
    if isinstance(node, FieldRef):
        return f"{node.name}_{node.idx}"

    # ── Time reference (integration only) ────────────────────────────────────
    if isinstance(node, TimeRef):
        if node.offset == 0:   n = node.name
        elif node.offset == -1: n = node.name + "_prev"
        else:                   n = node.name
        return f"vars.{n}[ i ]"

    # ── Dot product: ( a, b ) ────────────────────────────────────────────────
    if isinstance(node, Dot):
        return f"( {emit(node.a)}, {emit(node.b)} )"

    # ── Function call ────────────────────────────────────────────────────────
    if isinstance(node, Fn):
        args = ", ".join(emit(a) for a in node.args)
        return f"{node.name}( {args} )"

    # ── Binary op ────────────────────────────────────────────────────────────
    if isinstance(node, _B):
        mp = _PREC.get(node.op,5)
        l  = emit(node.left,  mp)
        r  = emit(node.right, mp)
        if node.op in ("-","/") and _p(node.right) <= mp:
            r = f"( {r} )"
        res = f"{l} {node.op} {r}"
        return f"( {res} )" if mp < pp else res

    # ── Unary op ─────────────────────────────────────────────────────────────
    if isinstance(node, _U):
        res = f"{node.op}{emit(node.operand, _PREC['u-'])}"
        return f"( {res} )" if _PREC["u-"] < pp else res

    # ── Comparison ───────────────────────────────────────────────────────────
    if isinstance(node, _C):
        mp = _PREC.get(node.op,3)
        res = f"{emit(node.left,mp)} {node.op} {emit(node.right,mp)}"
        return f"( {res} )" if mp < pp else res

    # ── Ternary ──────────────────────────────────────────────────────────────
    if isinstance(node, Ternary):
        res = f"( {emit(node.cond)} ) ? ( {emit(node.yes)} ) : ( {emit(node.no)} )"
        return f"( {res} )" if pp > 0 else res

    # ── _TimeRef2 (integration targets / rhs values) ─────────────────────────
    cls_name = type(node).__name__
    if cls_name == "_TimeRef2":
        return _emit_time_ref2(node)

    raise TypeError(f"emit: unknown node {type(node).__name__}")

def _emit_time_ref2(node) -> str:
    if node._offset == 0:
        arr = node._name
    elif node._offset == -1:
        arr = node._name + "_prev"
    else:
        arr = node._name
    return f"vars.{arr}[ i ]"
