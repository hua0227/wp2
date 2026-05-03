from __future__ import annotations

import unittest

import openmm_readtest_ligand_bound as readtest


class FakePlatform:
    names = ["CUDA", "OpenCL", "CPU"]

    def __init__(self, name: str) -> None:
        self._name = name

    def getName(self) -> str:
        return self._name

    @classmethod
    def getNumPlatforms(cls) -> int:
        return len(cls.names)

    @classmethod
    def getPlatform(cls, index: int) -> "FakePlatform":
        return cls(cls.names[index])

    @classmethod
    def getPlatformByName(cls, name: str) -> "FakePlatform":
        if name not in cls.names:
            raise Exception(f"No platform {name}")
        return cls(name)


class PlatformSelectionTests(unittest.TestCase):
    def test_choose_platform_prefers_opencl_when_available(self) -> None:
        platform = readtest.choose_platform(FakePlatform, preferred="OpenCL")

        self.assertEqual(platform.getName(), "OpenCL")


if __name__ == "__main__":
    unittest.main()
