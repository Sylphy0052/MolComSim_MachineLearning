import re
from enum import Enum
import linecache
import sys
import os
import numpy as np
from statistics import mean, variance, stdev, median
import math
from collections import deque

DNA_DICT = {
    11.31: {"diameter": 0.01, "dc": 63.944},
    9.19: {"diameter": 0.02, "dc": 42.187},
    7.46: {"diameter": 0.02, "dc": 27.833},
    6.06: {"diameter": 0.04, "dc": 18.363},
    4.92: {"diameter": 0.06, "dc": 12.115},
    4.00: {"diameter": 0.09, "dc": 7.993},
    3.25: {"diameter": 0.13, "dc": 5.273},
    2.64: {"diameter": 0.20, "dc": 3.479},
    2.14: {"diameter": 0.30, "dc": 2.295},
    1.18: {"diameter": 0.98, "dc": 0.695},
    0.67: {"diameter": 3.00, "dc": 0.227},
    2.20: {"diameter": 0.28, "dc": 2.416},
    1.42: {"diameter": 0.67, "dc": 1.013},
    10.58: {"diameter": 0.01, "dc": 55.931},
    1.00: {"diameter": 1.00, "dc": 0.500},
}

class Dat:
    def __init__(self, config_file_name):
        config_file_name = config_file_name
        self.config_header = []
        self.config_dict = {}

        self.initDict()
        self.parseConfigFile(config_file_name)

        self.input_data = InputData(self, self.config_dict)
        self.output_data = OutputData(self)

    def initDict(self):
        self.config_dict["moleculeParams"] = []
        self.config_dict["microtubuleParams"] = []

    def parseConfigFile(self, config_file_name):
        with open(config_file_name, 'r') as f:
            for line in f:
                if line[0] == '*' or line[0] == '\n':
                    continue
                line = line.rstrip()
                key, val = line.split(" ", 1)
                if not key in self.config_header:
                    self.config_header.append(key)

                if key in ["transmitter", "receiver"]:
                    self.config_dict[key] = NanoMachine(val)
                elif key in "intermediateNode":
                    self.config_dict[key] = IntermediateNode(val)
                elif key in "moleculeParams":
                    self.config_dict[key].append(MoleculeParams(val))
                elif key in "microtubuleParams":
                    self.config_dict[key].append(MicrotubuleParams(val))
                elif key in ["probDRail", "stepLengthX", "stepLengthY", "stepLengthZ"]:
                    self.config_dict[key] = float(val)
                elif key in "outputFile":
                    self.config_dict[key] = val
                elif key in ["useCollisions", "useAcknowledgements", "decomposing"]:
                    ret = False
                    if val ==str(1):
                        ret = True
                    self.config_dict[key] = ret
                else:
                    self.config_dict[key] = int(val)

    def getAllInfo(self):
        return [
                self.config_dict["mediumDimensionX"],
                self.output_data.analytical_model.r,
                self.config_dict["stepLengthX"],
                self.config_dict["moleculeParams"][0].num_of_molecules,
                self.config_dict["moleculeParams"][0].type_of_movement.value,
                # self.input_data.steps,
                ]

    def getColumn(self):
        return [
                "medium",
                "distance",
                "stepLength",
                "duplication",
                "movementType",
                "step"
                ]

class InputData:
    def __init__(self, dat, config_dict):
        input_file_name = "./result/batch_" + config_dict["outputFile"]
        self.steps = []
        self.ptime = []
        self.coll = []

        self.checkConfig(dat, input_file_name)
        self.parseData(dat, input_file_name)

        if dat.is_ptime:
            self.calcTxRxMean()
        if dat.is_coll:
            self.calcColl()

    def checkConfig(self, dat, input_file_name):
        length = len(linecache.getline(input_file_name, 1).rstrip().split(','))
        dat.is_ptime = False
        dat.is_coll = False

        if length == 5:
            dat.is_ptime = True
        elif length == 6:
            dat.is_coll = True
        elif length == 10:
            dat.is_ptime = True
            dat.is_coll = True

    def parseData(self, dat, input_file_name):
        inputdata = []
        with open(input_file_name, 'r') as f:
            for line in f:
                if ',' in line:
                    lines = line.rstrip().split(',')
                    inputdata.append(list(map(int, lines)))
                else:
                    inputdata.append(int(line))

        q = deque(zip(*inputdata))

        self.steps = np.sort(q.popleft())

        if dat.is_ptime:
            self.ptime.append(q.popleft())
            self.ptime.append(q.popleft())
            self.ptime.append(q.popleft())
            self.ptime.append(q.popleft())

        if dat.is_coll:
            self.coll.append(q.popleft())
            self.coll.append(q.popleft())
            self.coll.append(q.popleft())
            self.coll.append(q.popleft())
            self.coll.append(q.popleft())

    def calcTxRxMean(self):
        info_sum_time = sum(self.ptime[0])
        info_sum_num = sum(self.ptime[1])
        ack_sum_time = sum(self.ptime[2])
        ack_sum_num = sum(self.ptime[3])
        self.txrx_mean = float(info_sum_time) / info_sum_num + float(ack_sum_time) / ack_sum_num

    def calcColl(self):
        for i in range(len(self.coll)):
            self.coll[i] = sum(self.coll[i])

class OutputData:
    def __init__(self, dat):
        steps = dat.input_data.steps
        self.count = len(steps)
        self.minimum = int(steps[0])
        self.maximum = int(steps[-1])
        self.range_num = int(self.maximum / 1000) * 10
        if self.range_num == 0:
            self.range_num = int(self.maximum / 100)
        self.var = np.var(steps)
        self.std = np.std(steps)
        self.med = median(steps)
        self.mean = mean(steps)

        if dat.is_ptime:
            self.txrx_mean = dat.input_data.txrx_mean
        self.analytical_model = AnalyticalModel(dat.config_dict)

class AnalyticalModel:
    def __init__(self, config_dict):
        tx_position = config_dict["transmitter"].center_position
        rx_position = config_dict["receiver"].center_position
        step_length = config_dict["stepLengthX"]
        if step_length in DNA_DICT.keys():
            self.D = DNA_DICT[step_length]["dc"]
        else:
            self.D = 0.5

        self.r = self.calc_distance(tx_position, rx_position)
        self.L = int(config_dict["mediumDimensionX"] / 2)
        self.l = (int(config_dict["receiver"].size) * 2 - 1) / 2

        for molecule_param in config_dict["moleculeParams"]:
            # from transmitter to receiver
            if molecule_param.type_of_molecule == MoleculeType["INFO"]:
                if molecule_param.type_of_movement == MovementType["PASSIVE"]:
                    info_winf = self.calc_passive_rtt()
                else:
                    passive_winf = self.calc_passive_rtt()
                    info_winf = self.calc_active_rtt(passive_winf)

            # from receiver to transmitter
            elif molecule_param.type_of_molecule == MoleculeType["ACK"]:
                if molecule_param.type_of_movement == MovementType["PASSIVE"]:
                    ack_winf = self.calc_passive_rtt()
                else:
                    passive_winf = self.calc_passive_rtt()
                    ack_winf = self.calc_active_rtt(passive_winf)

        self.rtt = info_winf + ack_winf

    def calc_distance(self, tx, rx):
        x = (tx.x - rx.x) ** 2
        y = (tx.y - rx.y) ** 2
        z = (tx.z - rx.z) ** 2
        return math.sqrt(x + y + z)

    def calc_passive_rtt(self):
        winf = (self.r - self.l) * (2 * self.L ** 3 - self.l * self.r ** 2 - self.l ** 2 * self.r) / (2 * self.D * self.l * self.r)
        return winf

    def calc_active_rtt(self, passive_winf):
        V1 = 2 * 2 * self.r
        V = 4.0 / 3.0 * math.pi * self.L ** 3
        p = V1 / V
        va = 1
        vp = self.r / passive_winf
        ve = p * va + (1.0 - p) * vp
        winf = self.r / ve
        return winf

# 3次元位置情報
class Position:
    def __init__(self, args):
        self.x = int(args[0])
        self.y = int(args[1])
        self.z = int(args[2])

# transmitterとreceiverのクラス
class NanoMachine:
    def __init__(self, val):
        args = self.parse_val(val)
        self.center_position = Position(args[0:3])
        self.size = int(args[3])
        self.release_position = Position(args[4:7])

    def parse_val(self, val):
        val = [i for i in re.split(r"[,( )]", val) if i != '']
        return val

# 中間ノードクラス
class IntermediateNode:
    def __init__(self, val):
        args = self.parse_val(val)
        self.center_position = Position(args[0:3])
        self.size = int(args[3])
        self.info_release_position = Position(args[4:7])
        self.ack_release_position = Position(args[7:10])

    def parse_val(self, val):
        val = [i for i in re.split(r"[,( )]", val) if i != '']
        return val

# 分子の情報のクラス
class MoleculeParams:
    def __init__(self, val):
        args = self.parse_val(val)
        self.num_of_molecules = int(args[0])
        self.type_of_molecule = MoleculeType[args[1]]
        if self.type_of_molecule != MoleculeType["NOISE"]:
            self.type_of_movement = MovementType[args[2]]
            self.adaptive_change_number = int(args[3])
            try:
                self.size = float(args[4])
            except:
                pass
        else:
            try:
                self.size = float(args[2])
            except:
                pass

    def parse_val(self, val):
        val = [i for i in re.split(r"[ ]", val) if i != '']
        return val

# マイクロチューブのクラス
class MicrotubuleParams:
    def __init__(self, val):
        args = self.parse_val(val)
        self.start_position = Position(args[0:3])
        self.end_position = Position(args[3:6])

    def parse_val(self, val):
        val = [i for i in re.split(r"[,( )]", val) if i != '']
        return val

class MoleculeType(Enum):
    INFO = 0
    ACK = 1
    NOISE = 2

class MovementType(Enum):
    PASSIVE = 0
    ACTIVE = 1
