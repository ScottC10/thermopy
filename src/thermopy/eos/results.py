from dataclasses import dataclass

@dataclass(frozen=True)
class EoSResult:
    T : float
    P : float
    V : float
    Z : float
    phase : str | None
    eos_name : str
