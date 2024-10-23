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

#ifndef TESTRUNNINGVERTEXBASEDSIMULATIONSTUTORIAL_AB1_HPP_
#define TESTRUNNINGVERTEXBASEDSIMULATIONSTUTORIAL_AB1_HPP_

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

#include "StraightVertexMeshGenerator_ab1.hpp"

#include "StraightVertexMeshGenerator_ab2.hpp"		//added 9-17-21

#include "StraightVertexMeshGenerator_ab3.hpp"		//added 3-21-22

#include "CylindricalHoneycombVertexMeshGenerator.hpp"
/* The next header file defines a vertex-based {{{CellPopulation}}} class.*/
#include "VertexBasedCellPopulation.hpp"
/* The next header file defines a force law for describing the mechanical interactions
 * between neighbouring cells in the cell population, subject to each vertex.
 */

#include "CollagenCellLinkForce_ab.hpp"		 //added 7-14-21
#include "CollagenHookeForce_ab.hpp" 		//added 7-14-21
#include "SubstrateCollagenLinkForce_ab.hpp"		 //added 7-14-21
#include "IntraCellForce_ab.hpp"		 //added 7-14-21
#include "NagaiHondaForce.hpp"		 //added 7-14-21
#include "NagaiHondaForce_ab.hpp" //added 6-23-20
#include "SubstrateCellLinkForce_ab.hpp" //added 7-17-20
#include "SubstrateHookeForce_ab.hpp" 	//added 8-3-20

#include "CellProtrusiveContractileForce_ab.hpp"		 //added 9-17-21

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
class TestRunningVertexBasedSimulationsTutorial_ab1 : public AbstractCellBasedTestSuite
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
	bool durotaxis = true;		//added 8-11-23			

	unsigned numsubnode = 3;	//added 3-22-22		number of division between 0.5 lenth of substrate/collagen mesh	

	double cell_rest_length = 0.5;		//added 8-9-21
	//double collagen_rest_length = 0.5/3;		//added 8-9-21
	double collagen_rest_length = 0.5/numsubnode;		//added 3-22-22
	//double substrate_rest_length = 0.5/3;		//added 8-9-21
	double substrate_rest_length = 0.5/numsubnode;		//added 3-22-22



	double cell_spring_constant  = 50;		//added 8-9-21
	double collagen_spring_constant  = 1500;		//added 8-9-21
	double substrate_spring_constant  = 500;		//added 8-9-21
	double substrate_durotaxis_min_spring_constant = 350;	//added 8-11-23
	double substrate_durotaxis_max_spring_constant = 1000;		//added 8-11-23

	double collagen_cell_link_spring_constant  = 2500;		//added 8-9-21
	double substrate_collagen_link_spring_constant  = 1000;		//added 8-9-21

	
	double collagen_damping_constant  = 2;		//added 8-9-21
	double substrate_damping_constant  = 2;		//added 8-9-21
	bool variable_damping_constant = true;		//added 12-22-21
	double constant_cell_damping_constant  = 5;		//added 8-9-21, modified 12-22-21
	double min_cell_damping_constant = 5;		//added 12-22-21
	//double stdev_damping_constant = 1;		//added 12-22-21
	//double cell_damping_factor = 150;		//added 12-22-21


	double protrusive_force  = 6750;		//added 9-20-21

	

	double simulation_time = 20;		//added 11-24-21, 2 is 1000 timesteps fyi

	unsigned num_time_steps = 10000;		//added 12-10-21 max is equal to (simulation_time*500)

	double time_step = 0.0002;		//added 12-14-21 defaut is 0.002 hr for vertex based simulation


	bool threshold_based_division = true; 		//added 11-30-21 

	bool variable_threshold_length = true;		//added 12-20-21 if variable_threshold_length is true, then there will be the maximum threshold length for the cell mesh following a normal distribution. If stdev_div_length is 1, min will be max - 0.4. If std=2, min is max - 0.2

	double max_threshold_division_length = 0.675;		//added 12-20-21

	double constant_threshold_division_length = 0.95;	//added 9-22-21 ,modified 12-20-21. To be used in case variable_threshold_division is set 'false'

	//double stdev_div_length = 10;		//added 12-20-21


	bool islongcollagen = false;		//added 12-9-21

	//bool islongcollagen = true;		//added 12-9-21

	unsigned short_fiber_connectivity = 2;		//added 12-25-21

	bool update_adhesion1 = true;		//added 12-10-21

	double num_nodes_removed = 1; 		//added 12-25-21

	
	//double mean_cell_spring_constant  = 50;		//added 11-25-21	
	double stdev_cell_spring_constant  = 10;		//added 11-25-21
	
	unsigned cell_nodes = 4;// 36;		//added 12-19-21 actual nodes of cell are (cell_nodes+1). here it reresent number of cells

	double ND_factor_kcell = 100;		//added 12-19-21 this factor multiplies the normal distribution of spring stiffness within the cell mesh. 

	double ND_factor_tdl = 2;

	double ND_factor_dc = 75;
		
	double ND_stdev = 10;

	double ND_stdev_tdl = 5;

	double uv = (cell_nodes/2)+1;
	double s = ND_stdev;
	//std::normal_distribution<double> distribution(a, 1.0);
	std::vector<double> vect;
	//vect.push_back(this->mNodes.size());
	for (unsigned node_index = 0; node_index < cell_nodes+1; node_index++)
    	{
		double k_s = 1*(1/(s*sqrt(2*M_PI)))*exp(-(pow(((node_index+1-uv)/s),2.0))/2);
		//std::cout<<node_index+1-uv<<std::endl;
		//Node<SPACE_DIM>* p_node = this->GetNode(node_index);
		//p_node->AddElementStiffnessContribution(k_s);
		vect.push_back(k_s);
        	//std::cout<<p_node->rGetElementStiffnessContribution()<<std::endl;
    	}
	double max_ND = *max_element(vect.begin(), vect.end());
	double min_ND = *min_element(vect.begin(), vect.end());	
	
	double uv1 = (cell_nodes/2)+1;
	double s1 = ND_stdev_tdl;
	//std::normal_distribution<double> distribution(a, 1.0);
	std::vector<double> vect1;
	//vect.push_back(this->mNodes.size());
	for (unsigned node_index1 = 0; node_index1 < cell_nodes+1; node_index1++)
    	{
		double k_s1 = 1*(1/(s*sqrt(2*M_PI)))*exp(-(pow(((node_index1+1-uv1)/s1),2.0))/2);
		//std::cout<<node_index+1-uv<<std::endl;
		//Node<SPACE_DIM>* p_node = this->GetNode(node_index);
		//p_node->AddElementStiffnessContribution(k_s);
		vect1.push_back(k_s1);
        	//std::cout<<p_node->rGetElementStiffnessContribution()<<std::endl;
    	}
	double max_ND1 = *max_element(vect1.begin(), vect1.end());
	double min_ND1 = *min_element(vect1.begin(), vect1.end());

//********************assert(numElementsAcross % 4 == 0)*******************************************************************************************************		

	//StraightVertexMeshGenerator_ab1 gen1(8,1,4,false,false);    //added 7-1-21   1-D cell mesh
	
	//StraightVertexMeshGenerator_ab1 gen3(12,1,2,false,true,true,1);    //added 7-8-21   1-D collagen mesh. format - is substrate, iscollagen, islongcollagen, numnodesremoved 

	//StraightVertexMeshGenerator_ab1 gen2(16,1,0,true,false);    //added 7-1-21   1-D susbtrate mesh

	StraightVertexMeshGenerator_ab2 gen1(cell_nodes,1,4,false,false);    //added 7-1-21   1-D cell mesh
	
	//StraightVertexMeshGenerator_ab2 gen3(80,1,2,false,true,islongcollagen,num_nodes_removed,short_fiber_connectivity);    //added 7-8-21   1-D collagen mesh. format - is substrate, iscollagen, islongcollagen, numnodesremoved, shortfiberconnectivity  //commented 3-22-22 

	StraightVertexMeshGenerator_ab3 gen3(48,1,2,false,true,islongcollagen,num_nodes_removed,short_fiber_connectivity,numsubnode);    //added 7-8-21

	//StraightVertexMeshGenerator_ab2 gen2(8,1,0,true,false);    //added 7-1-21   1-D susbtrate mesh

	StraightVertexMeshGenerator_ab3 gen2(60,1,0,true,false,false,0,0,numsubnode);    //added 3-21-22   1-D susbtrate mesh
//**************************************************************************************************************************************************************
	//exit(0);

	//StraightVertexMeshGenerator_ab1 gen5(8,1,2,true,true,3);
	//exit(0);
	//exit(0);

	MutableVertexMesh<2,2>* p_mesh_ab1 = gen1.GetMesh();    //added 7-1-21   1-D cell mesh
	MutableVertexMesh<2,2>* p_mesh_ab2 = gen2.GetMesh();    //added 7-1-21   1-D substrate mesh
	MutableVertexMesh<2,2>* p_mesh_ab3 = gen3.GetMesh();    //added 7-8-21   1-D collagen mesh

	//exit(0);

        std::cout << "antman"<<std::endl;

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


	std::vector<CellPtr> cells1;    ////added 7-1-21

	//std::cout<<cells.size()<<std::endl;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

	CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator1;

	//cells_generator1.GenerateBasicRandom(cells1, p_mesh_ab1->GetNumElements(), p_transit_type);

	cells_generator1.GenerateBasicRandom_ab(constant_threshold_division_length, max_threshold_division_length,variable_threshold_length, ND_stdev_tdl, max_ND1, min_ND1, ND_factor_tdl, cells1,  p_mesh_ab1->GetNumElements(), p_transit_type);    //added 7-1-21, modified 12-20-21

	std::cout << "superman2\t" <<-stdev_cell_spring_constant<<std::endl;         
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

	VertexBasedCellPopulation<2> cell_population1(*p_mesh_ab1, cells1);    //added 7-1-21

	//std::cout << "superman333" << std::endl;
        /* We then pass the cell population into an {{{OffLatticeSimulation}}},
         * and set the output directory and end time. */

	//OffLatticeSimulation_ab<2> simulator_ab1(cell_population1, *p_mesh_ab2);		//added 7-1-21
	OffLatticeSimulation_ab1<2> simulator_ab2(cell_population1, *p_mesh_ab2, *p_mesh_ab3) ;		//added 7-8-21
	//exit(0);



	simulator_ab2.SetOutputDirectory("1-DVertexBasedMonolayer");//added 7-12-2
	//simulator_ab2.SetEndTime(2.0);		//added 7-12-21
	//simulator_ab2.SetEndTime(0.3);		//added 11-23-21	
	simulator_ab2.SetEndTime(simulation_time);		//added 7-12-21

        /*
         * For longer simulations, we may not want to output the results
         * every time step. In this case we can use the following method,
         * to print results every 50 time steps instead. As the default time step
         * used by the simulator (for vertex based simulations), is 0.02 hours, this method will cause the
         * simulator to print results every 6 minutes (i.e. 0.1 hours).
         */

	std::cout << "superman444" << std::endl;
        //simulator.SetSamplingTimestepMultiple(50);
	//simulator_ab.SetSamplingTimestepMultiple(50);//added 6-17-20
	//simulator_ab.SetSamplingTimestepMultiple(1);
	simulator_ab2.SetSamplingTimestepMultiple(1);		//added 7-12-21
	//exit(0);
	
	simulator_ab2.SetThresholdBasedCellDivision(threshold_based_division);		//added 11-30-21

	simulator_ab2.SetUpdateAdhesion1(update_adhesion1);		//added 12-10-21	

	simulator_ab2.SetNumTimeSteps(num_time_steps);		//added 12-10-21

	simulator_ab2.SetDt(time_step);		//added 12-14-21

	
        /* We must now create one or more force laws, which determine the mechanics of the vertices
        * of each cell in a cell population. For this test, we use one force law, based on the
        * Nagai-Honda mechanics, and pass it to the {{{OffLatticeSimulation}}}.
        * For a list of possible forces see subclasses of {{{AbstractForce}}}.
        * These can be found in the inheritance diagram, here, [class:AbstractForce AbstractForce].
        * Note that some of these forces are not compatible with vertex-based simulations see the specific class documentation for details,
        * if you try to use an incompatible class then you will receive a warning.
        */
        //MAKE_PTR(NagaiHondaForce<2>, p_force);

	MAKE_PTR(IntraCellForce_ab<2>, p_force_ab1); 	//added 7-12-21

	MAKE_PTR(CollagenCellLinkForce_ab<2>, p_force_ab2); 	//added 7-14-21	

	MAKE_PTR(CollagenHookeForce_ab<2>, p_force_ab3); 	//added 7-14-21

	MAKE_PTR(SubstrateCollagenLinkForce_ab<2>, p_force_ab4); 	//added 7-14-21
	
	MAKE_PTR(SubstrateHookeForce_ab<2>, p_force_ab5); //added 8-3-20	

	MAKE_PTR(SubstrateCellLinkForce_ab<2>, p_force_ab6); //added 7-17-20

	MAKE_PTR(CellProtrusiveContractileForce_ab<2>, p_force_ab7); 	//added 9-17-21
	


	p_mesh_ab1->SetRestLength(cell_rest_length);    //added 7-29-21   cell mesh spring rest length
	p_mesh_ab2->SetRestLength(substrate_rest_length);    //added 7-29-21   substrate mesh spring rest length
	p_mesh_ab3->SetRestLength(collagen_rest_length);    //added 7-29-21   collagen mesh spring rest length


	p_mesh_ab1->SetVariableSpringConstant(cell_nodes, stdev_cell_spring_constant, ND_factor_kcell);    //added 11-25-21   variable cell mesh spring constant 
	
	
	//exit(0);

	p_mesh_ab1->SetSpringConstant(cell_spring_constant);    //added 7-29-21   cell mesh spring constant
	
	if (durotaxis == false)		//added 8-12-23
	{
		p_mesh_ab2->SetSpringConstant(substrate_spring_constant);    //added 7-29-21   substrate mesh spring constant 
	}
	else
	{
		p_mesh_ab2->SetSpringConstantDurotaxis(substrate_spring_constant, substrate_durotaxis_min_spring_constant, substrate_durotaxis_max_spring_constant);    //added 8-11-23   substrate mesh spring constant for durotaxis
	}
	std::cout << "superman55" << std::endl;
	//exit(0);
	p_mesh_ab3->SetSpringConstant(collagen_spring_constant);    //added 7-29-21   collagen mesh spring constant	
	p_force_ab2->SetSpringConstant(collagen_cell_link_spring_constant);    //added 7-29-21   collagen-cell link spring constant
	p_force_ab4->SetSpringConstant(substrate_collagen_link_spring_constant);    //added 7-29-21   substrate-collagen link spring constant

	p_force_ab7->SetProtrusiveForce(protrusive_force);    //added 9-20-21   protrusive force on leader node

	//p_mesh_ab1->SetDampingConstant(cell_damping_constant);    //added 7-29-21   cell mesh Damping constant
	p_mesh_ab2->SetDampingConstant(substrate_damping_constant);    //added 7-29-21   substrate mesh Damping constant
	p_mesh_ab3->SetDampingConstant(collagen_damping_constant);    //added 7-29-21   collagen mesh Damping constant

	p_mesh_ab1->SetVariableDampingConstant(constant_cell_damping_constant, min_cell_damping_constant,variable_damping_constant, min_ND, max_ND, ND_stdev, ND_factor_dc);    //added 12-22-21   variable cell mesh damping constant 

	
	
			

	exit(0);		//uncommented 10-20-22 for plotting Kx, t_d_l, dc

	
	simulator_ab2.AddForce(p_force_ab1); 	//added 7-14-21  IntraCellForce
		
	simulator_ab2.AddForce(p_force_ab2); 	//added 7-14-21   CollagenCellLinkForce

	simulator_ab2.AddForce_ab1(p_force_ab3); 	//added 7-14-21   CollagenHookeForce

	simulator_ab2.AddForce_ab1(p_force_ab4); 	//added 7-14-21    SubstrateCollagenLinkForce   

	simulator_ab2.AddForce_ab(p_force_ab5); 	//added 7-14-21     SubstrateHookeForce  

	simulator_ab2.AddForce(p_force_ab7); 	//added 7-14-21  CellProtrusiveContractileForce

	

	//std::cout<<"simulator_ab2.GetDt()\t"<<simulator_ab2.GetDt()<<std::endl;
	//exit(0);

	std::cout << "superman555" << std::endl;
        /* A {{{NagaiHondaForce}}} assumes that each cell has a target area. The target areas of cells are used to determine pressure
         * forces on each vertex and eventually determine the size of each cell in the simulation. In order to assign target areas to cells
         * and update them in each time step we add a {{{SimpleTargetAreaModifier}}} to the simulation, which inherits from
         *  {{{AbstractTargetAreaModifier}}}.
         */
        //MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        //simulator.AddSimulationModifier(p_growth_modifier);
	//simulator_ab.AddSimulationModifier(p_growth_modifier);//added 6-17-20
	std::cout << "superman6" << std::endl;
	//exit(0);
        /* To run the simulation, we call {{{Solve()}}}. */
        //simulator.Solve();
	//simulator_ab.Solve();
	simulator_ab2.Solve();
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

#endif /* TESTRUNNINGVERTEXBASEDSIMULATIONSTUTORIAL_AB1_HPP_ */
