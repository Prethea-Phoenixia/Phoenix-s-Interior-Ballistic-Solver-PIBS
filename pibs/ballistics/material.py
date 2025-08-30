class Material:
    def __init__(self, rho: float, yield_strength: float, desc: str = ""):
        self.desc = desc
        self.yield_strength = yield_strength
        self.rho = rho
