"""
sph_codegen/dsl.py
===================
Symbolic expression tree nodes and statement records.

Every object that appears in a scheme expression is an Expr subclass.
Arithmetic operators return new Expr nodes — nothing is evaluated.
The code generator walks these trees later.

No monkey-patching. No execution side-effects. Pure data.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Union

# -- Expr base ---------------------------------------------------------------

class Expr:
    def __add__(self, o):  return _B("+",  self, _w(o))
    def __radd__(self, o): return _B("+",  _w(o), self)
    def __sub__(self, o):  return _B("-",  self, _w(o))
    def __rsub__(self, o): return _B("-",  _w(o), self)
    def __mul__(self, o):  return _B("*",  self, _w(o))
    def __rmul__(self, o): return _B("*",  _w(o), self)
    def __truediv__(self, o):  return _B("/", self, _w(o))
    def __rtruediv__(self, o): return _B("/", _w(o), self)
    def __neg__(self):    return _U("-", self)
    def __lt__(self, o):  return _C("<",  self, _w(o))
    def __le__(self, o):  return _C("<=", self, _w(o))
    def __gt__(self, o):  return _C(">",  self, _w(o))
    def __ge__(self, o):  return _C(">=", self, _w(o))

def _w(v):
    return v if isinstance(v, Expr) else Lit(float(v) if isinstance(v, int) else v)

@dataclass
class Lit(Expr):
    value: Union[int, float]

@dataclass
class _B(Expr):
    op: str; left: Expr; right: Expr

@dataclass
class _U(Expr):
    op: str; operand: Expr

@dataclass
class _C(Expr):
    op: str; left: Expr; right: Expr

@dataclass
class Fn(Expr):
    name: str; args: list

@dataclass
class Dot(Expr):
    a: Expr; b: Expr

@dataclass
class Ternary(Expr):
    cond: Expr; yes: Expr; no: Expr

@dataclass
class ConstRef(Expr):
    name: str

@dataclass
class FieldRef(Expr):
    name: str
    idx: str
    pset: str

@dataclass
class LocalRef(Expr):
    name: str

@dataclass
class TimeRef(Expr):
    name: str
    pset: str
    offset: int

# -- Statement records -------------------------------------------------------

@dataclass
class Stmt:
    kind: str           # "set" or "add"
    field: str
    idx: str
    pset: str
    expr: Expr
    condition: object = None

@dataclass
class LocalDecl:
    name: str
    expr: Expr
