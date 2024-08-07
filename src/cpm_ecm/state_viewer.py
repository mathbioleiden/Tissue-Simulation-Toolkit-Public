"""Model state viewer

Based on artist.py in TST-MD
"""
import logging
from pathlib import Path
from typing import List

from libmuscle import Instance
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from ymmsl import Operator

from tissue_simulation_toolkit.ecm.ecm import ParticleType
from tissue_simulation_toolkit.cpm_ecm.state_plotter import StatePlotter


_logger = logging.getLogger(__name__)


def main() -> None:
    logging.basicConfig(level=logging.DEBUG)
    logging.getLogger('PIL').setLevel(logging.INFO)
    logging.getLogger('matplotlib').setLevel(logging.INFO)

    instance = Instance({
        Operator.S: ['cpm_state_in', 'ecm_state_in']})

    while instance.reuse_instance():
        Lx = instance.get_setting('Lx', 'float')
        Ly = instance.get_setting('Ly', 'float')
        image_height = instance.get_setting('image_height', 'int')
        plotter = StatePlotter(Lx, Ly, image_height)

        state_output_interval = instance.get_setting('state_output_interval', 'int')
        mcs = instance.get_setting('mcs', 'int')

        for i in range(mcs):
            if i % state_output_interval == 0:
                cpm_state_msg = instance.receive('cpm_state_in')
                ecm_state_msg = instance.receive('ecm_state_in')

                particles = ecm_state_msg.data['particles']

                plotter.draw(
                        i,
                        particles['positions'].array,
                        particles['types'].array,
                        ecm_state_msg.data['bonds']['groups'].array,
                        cpm_state_msg.data['pde'].array,
                        cpm_state_msg.data['cpm'].array, save=False)

if __name__ == '__main__':
    main()
