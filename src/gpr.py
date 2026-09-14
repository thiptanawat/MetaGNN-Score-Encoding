"""Strict GPR parser: AND propagates missing evidence, OR retains measured alternatives."""
import re
import numpy as np
TOK = re.compile(r"\(|\)|[^\s()]+")

def evaluate(rule, gval, policy="conservative"):
    if policy not in {"conservative", "permissive"}:
        raise ValueError("Unknown GPR policy")
    tokens = TOK.findall(rule or "")
    if not tokens: return None
    pos = 0
    def atom():
        nonlocal pos
        if pos >= len(tokens): raise ValueError("Incomplete GPR expression")
        token = tokens[pos]; pos += 1
        if token == "(":
            value = or_expr()
            if pos >= len(tokens) or tokens[pos] != ")": raise ValueError("Unbalanced GPR parentheses")
            pos += 1
            return value
        if token in {"and", "or", ")"}: raise ValueError("Expected GPR gene operand")
        value = gval.get(token)
        return None if value is None else np.asarray(value, dtype=float)
    def and_expr():
        nonlocal pos
        value = atom()
        while pos < len(tokens) and tokens[pos] == "and":
            pos += 1; other = atom()
            if value is None or other is None:
                value = None if policy == "conservative" else (other if value is None else value)
            else:
                value = np.minimum(value, other) if policy == "conservative" else np.fmin(value, other)
        return value
    def or_expr():
        nonlocal pos
        value = and_expr()
        while pos < len(tokens) and tokens[pos] == "or":
            pos += 1; other = and_expr()
            if value is None: value = other
            elif other is not None: value = np.fmax(value, other)
        return value
    result = or_expr()
    if pos != len(tokens): raise ValueError("Unconsumed GPR tokens")
    return result
