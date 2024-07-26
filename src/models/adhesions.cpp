/*

Copyright 1996-2006 Roeland Merks
Copyright 2023 Netherlands eScience Center

This file is part of Tissue Simulation Toolkit.

Tissue Simulation Toolkit is free software; you can redistribute
it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

Tissue Simulation Toolkit is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Tissue Simulation Toolkit; if not, write to the Free
Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
02110-1301 USA

*/

#include <stdio.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include "adhesion_creation.hpp"
#include "cell.hpp"
#include "cpm_ecm/io.hpp"
#include "dish.hpp"
#include "graph.hpp"
#include "info.hpp"
#include "parameter.hpp"
#include "plotter.hpp"
#include "profiler.hpp"
#include "random.hpp"
#include "util/muscle3/settings.hpp"
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <libmuscle/libmuscle.hpp>
#include <ymmsl/ymmsl.hpp>

#include <iostream>
#include <vector>

using namespace std;

using libmuscle::Data;
using libmuscle::DataConstRef;
using libmuscle::Instance;
using libmuscle::Message;
using libmuscle::PortsDescription;
using libmuscle::StorageOrder;
using ymmsl::Operator;

extern Parameter par;

std::unique_ptr<Instance> instance;

INIT {
  try {
    // Define initial distribution of cells
    CPM->GrowInCells(par.n_init_cells, par.size_init_cells, par.subfield);
    CPM->ConstructInitCells(*this);

    // If we have only one big cell and divide it a few times
    // we start with a nice initial clump of cells.
    //
    // The behavior can be changed in the parameter file.
    for (int i = 0; i < par.divisions; i++) {
      CPM->DivideCells();
    }

    CPM->InitialiseEdgeList();
  } catch (const char *error) {
    cerr << "Caught exception\n";
    std::cerr << error << "\n";
    exit(1);
  }
}

TIMESTEP {
  try {
    static int i = 0;
    static Dish *dish = new Dish();
    static Info *info = new Info(*dish, *this);
    static Plotter plotter = Plotter(dish, this);

    if (i == 0) {
      // request creation of initial adhesions
      CellECMInteractions interactions;
      auto adh_zone = adhesion_zone(*(dish->CPM));
      interactions.change_type_in_area.change_area = adh_zone;
      interactions.change_type_in_area.num_particles =
          par.num_initial_adhesions;
      interactions.change_type_in_area.from_type = ParticleType::free;
      interactions.change_type_in_area.to_type = ParticleType::adhesion;
      auto data_mem = encode_cell_ecm_interactions(interactions);
      instance->send("cell_ecm_interactions_out", Message(i, data_mem.first));
    } else {
      // get any adhesion particle movements from CPM and send them out
      auto interactions = dish->CPM->GetCellECMInteractions();
      auto data_mem = encode_cell_ecm_interactions(interactions);
      instance->send("cell_ecm_interactions_out", Message(i, data_mem.first));
    }

    dish->CPM->ResetCellECMInteractions();

    if (i >= par.relaxation) {
      if (par.useopencl) {
        PROFILE(opencl_diff,
                dish->PDEfield->SecreteAndDiffuseCL(dish->CPM, par.pde_its);)
      } else {
        for (int r = 0; r < par.pde_its; r++) {
          dish->PDEfield->Secrete(dish->CPM);
          dish->PDEfield->Diffuse(1);
        }
      }
    }

    auto ecm_boundary_state_msg = instance->receive("ecm_boundary_state_in");
    auto ecm_boundary_state =
        decode_ecm_boundary_state(ecm_boundary_state_msg.data());
    dish->CPM->SetECMBoundaryState(ecm_boundary_state);

    PROFILE(amoebamove, dish->CPM->AmoebaeMove(dish->PDEfield);)

    if (instance->is_connected("state_out")) {
      if (i % instance->get_setting_as<int64_t>("state_output_interval") == 0) {
        std::cerr << "i = " << i << ", sending on state_out" << std::endl;
        auto *cpm_sigma = dish->CPM->getSigma()[0];
        Data cpm_state =
            Data::grid(cpm_sigma,
                       {static_cast<std::size_t>(dish->CPM->SizeX()),
                        static_cast<std::size_t>(dish->CPM->SizeY())},
                       {"x", "y"}, StorageOrder::last_adjacent);
        // auto const &pde = dish->PDEfield;
        // auto *pde_sigma = pde->getSigma()[0][0];
        // Data pde_state =
        //     Data::grid(pde_sigma,
        //                {static_cast<std::size_t>(pde->Layers()),
        //                 static_cast<std::size_t>(pde->SizeX()),
        //                 static_cast<std::size_t>(pde->SizeY())},
        //                {"layer", "x", "y"}, StorageOrder::first_adjacent);
        Data state = Data::dict("cpm", cpm_state); //Data state = Data::dict("cpm", cpm_state, "pde", pde_state);
        instance->send("state_out", Message(i, state));
      }
    }

    if (par.graphics && !(i % par.storage_stride)) {
      PROFILE(all_plots, plotter.Plot();)
      char title[400];
      snprintf(title, 399, "CellularPotts: %.2f hr",
               dish->PDEfield->TheTime() / 3600);
      info->Menu();
    }

    i++;
  } catch (const char *error) {
    cerr << "Caught exception\n";
    std::cerr << error << "\n";
    exit(1);
  } catch (std::exception const &e) {
    // ensure we crash if there's a problem, Qt swallows exceptions
    std::terminate();
  }
  PROFILE_PRINT
}

void PDE::Secrete(CellularPotts *cpm) {
  const double dt = par.dt;
  for (int x = 0; x < sizex; x++) {
    for (int y = 0; y < sizey; y++) {
      // inside cells
      if (cpm->Sigma(x, y)) {
        std::cout << "if" << std::endl; //sigma[0][x][y] += par.secr_rate[0] * dt;
      } else {
        // outside cells
        std::cout << "else" << std::endl; //sigma[0][x][y] -= par.decay_rate[0] * dt * sigma[0][x][y];
      }
    }
  }
  PROFILE_PRINT
}

void PDE::DerivativesPDE(CellularPotts *cpm, PDEFIELD_TYPE *derivs, int x,
                         int y) {}

int PDE::MapColour(double val) {
  return (((int)((val / ((val) + 1.)) * 100)) % 100) + 155;
}

void Plotter::Plot() {
  graphics->BeginScene();
  graphics->ClearImage();

  plotPDEDensity();
  plotCPMCellTypes();
  plotCPMLines();
  plotPDEContourLines();

  graphics->EndScene();
}

int main(int argc, char *argv[]) {
  PortsDescription ports(
      {{Operator::O_I, {"cell_ecm_interactions_out", "state_out"}},
       {Operator::S, {"ecm_boundary_state_in"}}});
  instance = std::make_unique<Instance>(0, nullptr, ports);

  instance->reuse_instance();
  set_parameters_from_settings(*instance);
  Seed(par.rseed);

  try {
    start_graphics(argc, argv);
  } catch (const char *error) {
    std::cerr << error << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "An unknown exception was caught" << std::endl;
    return 1;
  }

  // This is a hack, the whole model is really supposed to be inside a while
  // loop guarded by this statement. The architecture here won't allow that
  // and fortunately we don't need to actually run more than once, but
  // the function still needs to be called here for MUSCLE3 to shut down
  // correctly.
  instance->reuse_instance();

  return 0;
}
