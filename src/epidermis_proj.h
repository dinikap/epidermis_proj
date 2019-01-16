// -----------------------------------------------------------------------------
//
// Copyright (C) The BioDynaMo Project.
// All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------
#ifndef EPIDERMIS_PROJ_H_
#define EPIDERMIS_PROJ_H_

#include "biodynamo.h"

namespace bdm {

//Define a custom cell MyCell

BDM_SIM_OBJECT(MyCell, Cell) {
  BDM_SIM_OBJECT_HEADER(MyCell, Cell, 1, can_divide_, cell_type_);

public:
  MyCellExt() {}
  explicit MyCellExt(const std::array<double, 3>& position) : Base(position) {}
  /*
  In the epidermis, stem cells and TA cells can divide
  3 ways of division occurs
  1. Self-proliferating (i.e. SC -> SC + SC)
  2. Symmetric (i.e. SC -> TA + TA)
  3. Asymmetic (i.e. SC -> TA + SC)

  Q: should I use TMother?

  Code from bdm demo - tumor_concept
  */
  template <typename TMother>
  MyCellExt(const CellDivisionEvent& event, TMother* mother) : Base(event, mother) {
    can_divide_[kIdx] = mother->can_divide_[mother->kIdx];
  }

  //Q: daughter can have a different state from mother
  template <typename TDaughter>
  void EventHandler(const CellDivisionEvent& event, TDaughter* daughter) {
    Base::EventHandler(event, daughter);
  }

  void SetCanDivide(bool d) { can_divide_[kIdx] = d; }
  bool GetCanDivide() const { return can_divide_[kIdx]; }

  void SetCellType(int t) { cell_type_[kIdx] = t; }
  int GetCellType() { return cell_type_[kIdx]; }

private:
  vec<bool> can_divide_;
  vec<int> cell_type_;
};
/*
  Created own function to create cells
  so all cells will start from a fixed bound and then migrate
  based on the substance concentration

  MODIFIED to work with latest version of BDM
  ref: model_initialiser.h, CreateCellsRandom
*/

//template <typename Function, typename TResourceManager = ResourceManager<>>
template <typename Function, typename TSimulation = Simulation<>>
static void MyCellCreator(double min, double max, int num_cells, Function cell_builder) {
  auto* sim = TSimulation::GetActive();
  auto* rm = sim->GetResourceManager();
  auto* random = sim->GetRandom();
  // Determine simulation object type which is returned by the cell_builder
  using FunctionReturnType = decltype(cell_builder({0, 0, 0}));

  auto container = rm->template Get<FunctionReturnType>();
  container->reserve(num_cells);

  // so cells will be created at random only on the x and y axis
  // z axis is used to move cells to final resting position
  for (int i = 0; i < num_cells; i++) {
    double x = random->Uniform(min, max);
    double y = random->Uniform(min, max);
    //stop cells from moving in the z axis when generated
    double z = 0;
    auto new_simulation_object = cell_builder({x, y, z});
    container->push_back(new_simulation_object);
  }
  container->Commit();
}

//First type of cell: Stem cells
//Stem cells will divide and grow to produce TA cells
struct StemCell : public BaseBiologyModule {
  StemCell() : BaseBiologyModule(gAllEventIds){}

  template <typename TEvent, typename TBm>
  StemCell(const TEvent& event, TBm* other, uint64_t newoid = 0) {}

  template <typename T, typename TSimulation = Simulation<>>
  void Run(T* cell) {
    //auto* sim = TSimulation::GetActive();
    // auto* rm = sim->GetResourceManager();
    // static auto* kDg = rm->GetDiffusionGrid(calcium);
    auto&& daughter = cell->Divide();
    if (cell->GetDiameter() < 5 && cell->GetCellType() == 1) {
      //self-proliferation -> divide to itself
      daughter->SetCellType(1);
      daughter->SetCanDivide(true);
    }else if (cell->GetDiameter() < 8 && cell->GetCellType() == 1) {
      //Symmetric division -> TA + TA
      daughter->SetCellType(2);
      daughter->SetCanDivide(true);
    }else
      cell->SetCanDivide(false); //inactivate cell
  }
  BDM_CLASS_DEF_NV(StemCell, 1);
};

//TA cells
struct TransitAmplifying : public BaseBiologyModule {
  TransitAmplifying() : BaseBiologyModule(gAllEventIds){}

  template <typename TEvent, typename TBm>
  TransitAmplifying(const TEvent& event, TBm* other, uint64_t newoid = 0) {}

  template <typename T, typename TSimulation = Simulation<>>
  void Run(T* cell) {
    //auto* sim = TSimulation::GetActive();
    // auto* rm = sim->GetResourceManager();
    // static auto* kDg = rm->GetDiffusionGrid(calcium);
    auto&& daughter = cell->Divide();
    if (cell->GetDiameter() < 8 && cell->GetCellType() == 2) {
      //self-proliferation -> divide to itself
      daughter->SetCellType(2);
      daughter->SetCanDivide(true);
    }else if (cell->GetDiameter() < 10 && cell->GetCellType() == 2) {
      //Symmetric division -> GA + GA
      daughter->SetCellType(3);
      daughter->SetCanDivide(true);
    }else
      cell->SetCanDivide(false); //inactivate cell
  }
  BDM_CLASS_DEF_NV(TransitAmplifying, 1);
};

//Differentiated cells
struct DifferentiatedCell : public BaseBiologyModule {
  DifferentiatedCell() : BaseBiologyModule(gAllEventIds){}

  template <typename TEvent, typename TBm>
  DifferentiatedCell(const TEvent& event, TBm* other, uint64_t new_oid = 0) {}

  template <typename T, typename TSimulation = Simulation<>>
  void Run(T* cell){
    //auto* sim = TSimulation::GetActive();
    // auto* rm = sim->GetResourceManager();
    // static auto* kDg = rm->GetDiffusionGrid(calcium);
    if (cell->GetDiameter() > 10){
      cell->SetCellType(3);
    }
  }

private:
  BDM_CLASS_DEF_NV(DifferentiatedCell, 1);
};

// Define compile time parameter
BDM_CTPARAM() {
  BDM_CTPARAM_HEADER();
  using SimObjectTypes = CTList<MyCell>;  // use MyCell object

  // Override default BiologyModules for Cell
  BDM_CTPARAM_FOR(bdm, MyCell) {
    using BiologyModules = CTList<StemCell, TransitAmplifying, DifferentiatedCell>;
  };
};

inline int Simulate(int argc, const char** argv) {
  //set parameters
  auto set_param = [](auto* param) {
    param->bound_space_ = true;
    param->min_bound_ = 0;
    param->max_bound_ = 250;  //TODO double check param bounds for epidermis
    param->run_mechanical_interactions_ = true; //or false?
  };

  Simulation<> simulation(argc, argv, set_param);

  // Define initial model - in this example: single cell at origin
  //auto* rm = simulation.GetResourceManager();
  //auto&& cell = rm->New<Cell>(30);

  auto* param = simulation.GetParam();
  auto* random = simulation.GetRandom();
  cout << "Random seed " << random << endl;
  auto construct_stem = [](const std::array<double, 3>& position) {
    MyCell cell(position);
    cell.SetDiameter(2);
    cell.AddBiologyModule(StemCell());
    cell.SetCellType(1);
    return cell;
  };
  cout << "Stem cells created" << endl;
  MyCellCreator(param->min_bound_, param->max_bound_, 200, construct_stem);

  // Run simulation for one timestep
  simulation.GetScheduler()->Simulate(1);

  std::cout << "Simulation completed successfully!" << std::endl;
  return 0;
}

}  // namespace bdm

#endif  // EPIDERMIS_PROJ_H_
