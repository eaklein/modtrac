# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 10:32:02 2021

@author: Avram
"""
import time
import math
import numpy as np
import pandas as pd
from pathlib import Path
import glob

from .utilities import t_to_d


class PtracMod:
    """A class to store ptrac information for moderator studies."""

    def __init__(self, key, thick, mode="n", folder=Path(".")):
        """Initialize TargetData class."""
        self.key = key
        self.mode = mode
        self.thick = thick
        self.data = pd.DataFrame()
        self.data_axial = pd.DataFrame()
        self.geometry = {}
        self.color = "black"
        self.label = self.key + "-" + self.thick
        self.nps = 0
        self.seeds = []
        self.folder = folder + self.key
        self.file_dat = (
            Path(self.folder) / self.thick / (
                self.key + "-" + self.thick + ".dat")
        )
        self.file_outp = ""
        self.file_ptrac = []
        self.file_hdf = ""

    def get_MCNP_params(self):
        """Read in run parameters for MCNP run."""
        print(f"\nGetting parameters for {self.label}")
        # self.read_params()
        self.read_outp()

    def read_params(self):
        """Read custom-written README file to get run parameters."""
        params = {
            "type_mod": str,
            "shape_mod": str,
            "r_mod": float,
            "z_mod": float,
            "th_mod": float,
            "type_shield": str,
            "shape_shield": str,
            "r_shield": float,
            "z_shield": float,
        }
        try:
            with open(self.folder / "params.txt") as f:
                next(f)
                for line in f:
                    if len(line.split()) == 2:
                        param, val = line.split()
                        self.geometry[param] = params[param](val.strip())
            print(f"Read params.txt file for key: {self.key}")
        except:
            print(f"Error in reading params.txt file for key: {self.key}")
        r_dt = 5.0
        self.geometry["area_mod"] = np.pi * (
            (self.geometry["r_mod"] + r_dt) ** 2 - (r_dt) ** 2
        )

    def read_outp(self):
        """Read nps and seed from corresponding output files."""
        files_outp = [
            str(x)
            for x in glob.glob(self.folder + "/" + self.thick + "/out*",
                               recursive=True)
        ]
        files_ptrac = []
        nps_total = 0
        if len(files_outp) == 1:
            print("Found 1 output file for key: "
                  f"{self.key} (mode: {self.mode})")
        elif len(files_outp) > 1:
            print(
                f"Found {len(files_outp)} output files for key: "
                "{self.key} (mode: {self.mode})"
            )
        else:
            print(
                "No output files found for key: "
                f"{self.key} (mode: {self.mode}) in {self.folder}"
            )
        # determine total number of particles simulated
        seeds = []
        for file in files_outp:
            nps_file = 0
            with open(file) as f:
                for line in f:
                    if line.strip().startswith("* Random Number Seed"):
                        seed = int(line.split("=")[1].split("*")[0].strip())
                        seeds.append(seed)
                        # print(f'Random Seed: {seed}')
                    elif line.strip().startswith("run terminated"):
                        if line.strip() == "run terminated by tty interrupt.":
                            pass
                        elif line.strip().endswith("computer time."):
                            pass
                        else:
                            nps_file = int(
                                line.split("when")[1].split(
                                    "particle")[0].strip()
                            )
                            nps_total += nps_file
                            # print(f'NPS: {nps:1.2e}')
                    if line.strip().startswith("binary file"):
                        file_ptrac = line.split("binary file ")[1].split(
                            " ")[0]
                        files_ptrac.append(file_ptrac)
                # raise error if outp file terminated prematurely
                # and request user input
                if nps_file == 0:
                    print(
                        f"Error encountered in {file}. "
                        "Run terminated prematurely."
                    )
                    nps_last = 0
                    # read outp
                    with open(file) as f:
                        for line in f:
                            if line.strip().startswith("dump no."):
                                nps_dump = int(
                                    line.split("nps =")[1].split(
                                        "coll")[0].strip()
                                )
                                if nps_dump > nps_last:
                                    nps_last = nps_dump
                    print(f"Assuming {nps_last:1.2e} particles simulated.")
                    nps_file = nps_last
                    # nps = click.prompt('Manually input number of particles '
                    #                    f'simulated for {self.label}.',
                    #                    type=int, default=nps_last)
                    nps_total += nps_file
        self.file_ptrac = files_ptrac
        self.file_outp = files_outp
        self.seeds = seeds
        self.nps = nps_total

    def read_ptrac(self, mode="n", l_tof=2.00, rad=10.0, save=""):
        """Read in data and calculate position at extended distance."""
        # try:
        #     filename = self.folder + self.key + '-' + self.mode + '.csv'
        #     print(f'Trying to read in {filename}...', end='\r')
        #     df = pd.read_csv(filename)
        #     print('Done!')
        # except:
        filename = Path(self.folder) / self.thick / (
            self.key + "-" + self.thick)
        print(f"Attempting to process data for {self.key}-{self.thick}...",
              end="\r")
        time.sleep(1)
        try:
            self.file_hdf = (
                Path(self.folder) / self.thick / (
                    self.key + "-" + self.thick + ".h5")
            )
            self.read_hdf()
            print("HDF5 file processed!")
        except:
            print("Unable to read HDF5 file. Processing manually...")
            # for ptrac in self.file_ptrac:
            #     print(str(filename) + '-' + ptrac + '.dat')
            # if not Path(filename).is_file():
            #     print(f'{filename} does not exist.')
            #     return
            # print(f'Reading {filename}')
            # df = pd.concat([pd.read_csv(
            #     str(filename) + '-' + ptrac + '.dat',
            #     header=None, sep=" ", on_bad_lines='skip'
            #     ) for ptrac in self.file_ptrac])
            df = pd.read_csv(
                str(filename) + ".dat", header=None, sep=" ",
                on_bad_lines="skip"
            )
            # ensure matches ptrac data output in track.py
            df.columns = ["t", "E", "x", "y", "z", "u", "v", "w",
                          "weight", "coll"]
            df.t /= 100  # convert shakes to us
            df.E *= 1e6  # convert MeV to eV
            t_offset = 0.0  # in us
            if mode == "n":
                # convert mod time to equivalent TOF distance; time (s), E (eV)
                df["d"] = 1000 * t_to_d((df.t - t_offset) / 1e6, df.E)
            # radial dist from TOF axis at source-moderator assembly
            df["r"] = np.sqrt(df.x**2 + df.y**2)
            df["phi"] = math.atan2(df.y, df.x)
            # calculate n position at extended distance
            df["rad2"] = (df.x + (l_tof / df.w) * df.u) ** 2 + (df.y + (l_tof / df.w) * df.v) ** 2
            self.data = df
            self.data_axial = df.loc[df.rad2 <= rad**2].copy(deep=True)
            if save == "csv":
                df.to_csv(
                    Path(self.folder) / (self.key + "-" + self.thick + ".csv")
                )  # write .csv to file
            elif save.upper() == "HDF":
                if not Path(self.file_hdf).is_file():
                    self.write_hdf()
                else:
                    print("HDF5 file already exists for "
                          f"{self.key}-{self.thick}")

    def calc_axial(self, rad=10.0):
        """Acquire data for TOF-axial particles at given distance."""
        df = self.data.copy()
        self.data_axial = df.loc[df.rad2 <= rad**2]

    def write_csv(self):
        """Store processed ptrac as csv."""
        self.data.to_csv(self.folder + self.key + "-" + self.mode + ".csv")

    def read_csv(self):
        """Read csv of processed ptrac."""
        self.data = pd.read_csv(self.folder + self.key + "-" +
                                self.mode + ".csv")

    def read_hdf(self):
        """Read hdf of processed ptrac."""
        self.data = pd.read_hdf(self.file_hdf, key="data")
        self.data_axial = pd.read_hdf(self.file_hdf, key="data_axial")

    def write_hdf(self):
        """Store processed ptrac as hdf."""
        store = pd.HDFStore(self.file_hdf)
        store.append("data", self.data)
        store.append("data_axial", self.data_axial)
        store.close()
