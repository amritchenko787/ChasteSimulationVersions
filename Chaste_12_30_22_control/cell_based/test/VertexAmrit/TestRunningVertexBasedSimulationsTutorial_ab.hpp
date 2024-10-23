/*

Copyright (c) 2005-2018, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/
/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */

#ifndef TESTRUNNINGVERTEXBASEDSIMULATIONSTUTORIAL_AB_HPP_
#define TESTRUNNINGVERTEXBASEDSIMULATIONSTUTORIAL_AB_HPP_

/*
 * = Examples showing how to create, run and visualize vertex-based simulations =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how Chaste can be used to create, run and visualize vertex-based simulations.
 * Full details of the mechanical model proposed by T. Nagai and H. Honda ("A dynamic cell model for
 * the formation of epithelial tissues", Philosophical Magazine Part B 81:699-719).
 *
 * EMPTYLINE
 *
 * == The test ==
 *
 * EMPTYLINE
 *
 * As in previous cell-based Chaste tutorials, we begin by including the necessary header files.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

/* The remaining header files define classes that will be used in the cell-based
 * simulation. We have encountered some of these header files in previous cell-based
 * Chaste tutorials. */
#include "CellsGenerator.hpp"
#include "OffLatticeSimulation_ab.hpp"
#include "OffLatticeSimulation_ab1.hpp"    //added 7-8-21
#include "OffLatticeSimulation.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
/* The next header file defines the cell cycle model. */
#include "UniformG1GenerationalCellCycleModel.hpp"
/* The next two header files define a helper class for generating suitable meshes: one planar and one periodic. */
#include "HoneycombVertexMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator_ab.hpp"
#include "StraightVertexMeshGenerator_ab.hpp"
#include "CylindricalHoneycombVertexMeshGenerator.hpp"
/* The next header file defines a vertex-based {{{CellPopulation}}} class.*/
#include "VertexBasedCellPopulation.hpp"
/* The next header file defines a force law for describing the mechanical interactions
 * between neighbouring cells in the cell population, subject to each vertex.
 */

#include "IntraCellForce_ab.hpp"
#include "NagaiHondaForce.hpp"
#include "NagaiHondaForce_ab.hpp" //added 6-23-20
#include "SubstrateCellLinkForce_ab.hpp" //added 7-17-20
#include "SubstrateHookeForce_ab.hpp" //added 8-3-20

/* This force law assumes that cells possess a "target area" property which determines the size of each
 * cell in the simulation. In order to assign target areas to cells and update them in each time step, we need
 * the next header file.
 */
#include "SimpleTargetAreaModifier.hpp"
/* The next header file defines a boundary condition for the cells.*/
#include "PlaneBoundaryCondition.hpp"
/* The next header file defines a cell killer, which specifies how cells are removed from the simulation.*/
#include "PlaneBasedCellKiller.hpp"

/* Finally, we include a header that enforces running this test only on one process. */
#include "FakePetscSetup.hpp"

/* Next, we define the test class, which inherits from {{{AbstractCellBasedTestSuite}}}
 * and defines some test methods.
 */
class TestRunningVertexBasedSimulationsTutorial_ab : public AbstractCellBasedTestSuite
{
public:
    /* EMPTYLINE
    *
    * == Test 1 - a basic vertex-based simulation ==
    *
    * EMPTYLINE
    *
    * In the first test, we run a simple vertex-based simulation, in which we create a monolayer
    * of cells, using a mutable vertex mesh. Each cell is assigned a stochastic cell-cycle model.
    */
    void TestMonolayer()
    {
        /* First, we generate a vertex mesh. To create a {{{MutableVertexMesh}}}, we can use
        * the {{{HoneycombVertexMeshGenerator}}}. This generates a honeycomb-shaped mesh,
        * in which all nodes are equidistant. Here the first and second arguments
        * define the size of the mesh - we have chosen a mesh that is 2 elements (i.e.
        * cells) wide, and 2 elements high.
        */
        HoneycombVertexMeshGenerator_ab gen(2, 2);
	HoneycombVertexMeshGenerator_ab gen4(2, 2);
	//StraightVertexMeshGenerator_ab gen1(2,1,1,false);    //added 7-1-21   1-D cell mesh
	//StraightVertexMeshGenerator_ab gen2(3,1,0,false);    //added 7-1-21   1-D susbtrate mesh
	//StraightVertexMeshGenerator_ab gen3(3,1,0.5,true);    //added 7-8-21   1-D collagen mesh, trial mesh just to see if it can be passed onto offlatticesimulation_ab1
	//exit(0);
	MutableVertexMesh<2,2>* p_mesh_ab = gen.GetMesh();
	MutableVertexMesh<2,2>* p_mesh_ab1 = gen4.GetMesh();    //added 7-1-21   1-D cell mesh
	//MutableVertexMesh<2,2>* p_mesh_ab2 = gen2.GetMesh();    //added 7-1-21   1-D substrate mesh
	//MutableVertexMesh<2,2>* p_mesh_ab3 = gen3.GetMesh();    //added 7-8-21   1-D collagen mesh
	p_mesh_ab->Scale(0.2,0.2);	//added 4-23-21
        std::cout << "antman"<<std::endl;
        HoneycombVertexMeshGenerator generator(2, 2);    // Parameters are: cells across, cells up
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
	//std::cout<<"star"<<std::endl;	
	//std::cout<<*p_mesh<<std::endl;
 	//std::cout<<p_mesh<<std::endl;
	//std::cout<<&p_mesh<<std::endl;
	//std::cout<<"max"<<std::endl;
        /* Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
        * To do this, we use the `CellsGenerator` helper class, which is templated over the type
        * of cell model required (here {{{UniformG1GenerationalCellCycleModel}}})
        * and the dimension. We create an empty vector of cells and pass this into the
        * method along with the mesh. The second argument represents the size of that the vector
        * {{{cells}}} should become - one cell for each element, the third argument specifies
        * the proliferative type of the cell. */
        std::vector<CellPtr> cells;   //@Cell.hpp- typedef boost::shared_ptr<Cell> CellPtr;

	std::vector<CellPtr> cells1;    ////added 7-1-21

	//std::cout<<cells.size()<<std::endl;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
	CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator1;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

	cells_generator1.GenerateBasicRandom(cells1, p_mesh_ab1->GetNumElements(), p_transit_type);    //added 7-1-21  
	std::cout << "superman2" << std::endl;         
	//exit(0);
        
        //for (std::vector<CellPtr>::iterator it = cells.begin() ; it != cells.end(); ++it)
        //{ std::cout << ' ' << *it;
        //std::cout << '\n';
        //std::cout <<"superman"<<std::endl;}
	//std::cout<<cells.size()<<std::endl;
        /* Now we have a mesh and a set of cells to go with it, we can create a {{{CellPopulation}}}.
        * In general, this class associates a collection of cells with a mesh.
        * For this test, because we have a {{{MutableVertexMesh}}}, we use a particular type of
        * cell population called a {{{VertexBasedCellPopulation}}}.
        */
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
	VertexBasedCellPopulation<2> cell_population1(*p_mesh_ab1, cells1);    //added 7-1-21
	//std::cout << "superman333" << std::endl;
        /* We then pass the cell population into an {{{OffLatticeSimulation}}},
         * and set the output directory and end time. */
        OffLatticeSimulation_ab<2> simulator_ab(cell_population, *p_mesh_ab); //added 6-11-20
	//OffLatticeSimulation_ab<2> simulator_ab1(cell_population1, *p_mesh_ab2);		//added 7-1-21
	//OffLatticeSimulation_ab1<2> simulator_ab2(cell_population1, *p_mesh_ab2, *p_mesh_ab3) ;		//added 7-8-21
	//exit(0);
	OffLatticeSimulation<2> simulator(cell_population);


	//simulator_ab2.SetOutputDirectory("1-DVertexBasedMonolayer");//added 7-12-21
        simulator_ab.SetOutputDirectory("VertexBasedMonolayer");//added 6-17-20
	simulator.SetOutputDirectory("VertexBasedMonolayer");

	//simulator_ab2.SetEndTime(1.0);		//added 7-12-21
        simulator_ab.SetEndTime(1.0);//added 6-17-20
	simulator.SetEndTime(1.0);

        /*
         * For longer simulations, we may not want to output the results
         * every time step. In this case we can use the following method,
         * to print results every 50 time steps instead. As the default time step
         * used by the simulator (for vertex based simulations), is 0.02 hours, this method will cause the
         * simulator to print results every 6 minutes (i.e. 0.1 hours).
         */

	std::cout << "superman444" << std::endl;
        simulator.SetSamplingTimestepMultiple(50);
	//simulator_ab.SetSamplingTimestepMultiple(50);//added 6-17-20
	simulator_ab.SetSamplingTimestepMultiple(1);
	//simulator_ab2.SetSamplingTimestepMultiple(1);		//added 7-12-21

        /* We must now create one or more force laws, which determine the mechanics of the vertices
        * of each cell in a cell population. For this test, we use one force law, based on the
        * Nagai-Honda mechanics, and pass it to the {{{OffLatticeSimulation}}}.
        * For a list of possible forces see subclasses of {{{AbstractForce}}}.
        * These can be found in the inheritance diagram, here, [class:AbstractForce AbstractForce].
        * Note that some of these forces are not compatible with vertex-based simulations see the specific class documentation for details,
        * if you try to use an incompatible class then you will receive a warning.
        */
        //MAKE_PTR(NagaiHondaForce<2>, p_force);

	//MAKE_PTR(IntraCellForce_ab<2>, p_force_ab3); 	//added 7-12-21
	MAKE_PTR(NagaiHondaForce_ab<2>, p_force_ab); //added 6-23-20
	MAKE_PTR(SubstrateCellLinkForce_ab<2>, p_force_ab1); //added 7-17-20
	MAKE_PTR(SubstrateHookeForce_ab<2>, p_force_ab2); //added 8-3-20


        //simulator.AddForce(p_force);
	//simulator_ab.AddForce(p_force);//added 6-17-20
	//simulator_ab2.AddForce(p_force_ab3); 	//added 7-12-21
	simulator_ab.AddForce(p_force_ab);//added 6-23-20
	simulator_ab.AddForce(p_force_ab1);//added 6-23-20
	simulator_ab.AddForce_ab(p_force_ab2);//added 8-3-20
	std::cout << "superman555" << std::endl;
        /* A {{{NagaiHondaForce}}} assumes that each cell has a target area. The target areas of cells are used to determine pressure
         * forces on each vertex and eventually determine the size of each cell in the simulation. In order to assign target areas to cells
         * and update them in each time step we add a {{{SimpleTargetAreaModifier}}} to the simulation, which inherits from
         *  {{{AbstractTargetAreaModifier}}}.
         */
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);
	simulator_ab.AddSimulationModifier(p_growth_modifier);//added 6-17-20
	std::cout << "superman6" << std::endl;
	//exit(0);
        /* To run the simulation, we call {{{Solve()}}}. */
        //simulator.Solve();
	simulator_ab.Solve();
	//simulator_ab2.Solve();
	std::cout << "superman7" << std::endl;
        /* The next two lines are for test purposes only and are not part of this tutorial. If different simulation input parameters are being explored
         * the lines should be removed.*/
        //TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 4u);
        //TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 1.0, 1e-10);
    }

    /*
    * EMPTYLINE
    *
    * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
    * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dVertexCells /tmp/$USER/testoutput/VertexBasedMonolayer/results_from_time_0}}}.
    * We may have to do: {{{javac Visualize2dVertexCells.java}}} beforehand to create the
    * java executable.
    *
    * EMPTYLINE
    *
    * == Test 2 - introducing periodicity, boundaries and cell killers ==
    *
    * EMPTYLINE
    *
    * In the second test, we run a simple vertex-based simulation, in which we create a monolayer
    * of cells in a periodic geometry, using a cylindrical vertex mesh. We also include a fixed
    * boundary which cells can't pass through and a cell killer which removes cells once they leave
    * a region. As before each cell is assigned a stochastic cell-cycle model.
    */


};

#endif /* TESTRUNNINGVERTEXBASEDSIMULATIONSTUTORIAL_AB_HPP_ */
