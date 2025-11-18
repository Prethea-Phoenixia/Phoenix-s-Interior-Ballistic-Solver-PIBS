from __future__ import annotations
import json
from . import JSONable


class Material(JSONable):
    def __init__(self, density: float, yield_strength: float, desc: str = ""):
        self.desc = desc
        self.yield_strength = yield_strength
        self.density = density

    def to_json(self) -> str:
        return json.dumps(
            {
                "desc": self.desc,
                "yield_strength": self.yield_strength,
                "density": self.density,
            },
            ensure_ascii=False,
        )

    @classmethod
    def from_json(cls, json_dict: dict) -> Material:
        return cls(**json_dict)
