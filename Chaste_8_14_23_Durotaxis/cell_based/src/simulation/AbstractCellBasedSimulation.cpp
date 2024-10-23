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

#include <cmath>
#include <iostream>
#include <fstream>
#include <set>

#include "AbstractCellBasedSimulation.hpp"
#include "CellBasedEventHandler.hpp"
#include "LogFile.hpp"
#include "ExecutableSupport.hpp"
#include "AbstractPdeModifier.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::AbstractCellBasedSimulation(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                              bool deleteCellPopulationInDestructor,
                                              bool initialiseCells)
    : mIsCol(true),		  //added 7-16-21
      mDt(DOUBLE_UNSET),
      mEndTime(DOUBLE_UNSET),  // hours - this is set later on
      mrCellPopulation(rCellPopulation),
      mDeleteCellPopulationInDestructor(deleteCellPopulationInDestructor),
      mInitialiseCells(initialiseCells),
      mNoBirth(false),
      mUpdateCellPopulation(true), //Modified on 7-16-21 from mUpdateCellPopulation(true)	
      mOutputDirectory(""),
      mSimulationOutputDirectory(mOutputDirectory),
      mNumBirths(0),
      mNumDeaths(0),
      mOutputDivisionLocations(false),
      mOutputCellVelocities(true),
      mSamplingTimestepMultiple(1)
{
    // Set a random seed of 0 if it wasn't specified earlier
    RandomNumberGenerator::Instance();

    if (mInitialiseCells)
    {
        mrCellPopulation.InitialiseCells();
    }

    /*
     * Specify a default time step to use, which may depend on the cell population type.
     * Note that the time step may be reset using SetDt().
     */
    //mDt = rCellPopulation.GetDefaultTimeStep();	commented 12-14-21
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::~AbstractCellBasedSimulation()
{
    if (mDeleteCellPopulationInDestructor)
    {
        delete &mrCellPopulation;
    }
}



//added 9-21-21
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::DoCellBirth_ab()
{
    	//std::cout<<"inside_docellbirth_ab"<<std::endl;
    if (mNoBirth)
    {
        return 0;
    }

    unsigned num_births_this_step = 0;

    // Iterate over all cells, seeing if each one can be divided
    for (typename AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
	//std::cout<<"inside_docellbirth_loop"<<std::endl;
	double cell_size = mrCellPopulation.GetSizeOfCell_ab(*cell_iter);

	bool isleader = mrCellPopulation.IsLeaderCell_ab(*cell_iter); 
	//old_cell_locations[*cell_iter] = mrCellPopulation.GetLocationOfCellCentre(*cell_iter);
	//VertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = mrCellPopulation.GetElementCorrespondingToCell_ab(*cell_iter);		//added 12-17-21
	//Node<SPACE_DIM>* p_node_A = p_element->GetNode(0);
	//Node<SPACE_DIM>* p_node_B = p_element->GetNode(1);
	//bool test_A = p_node_A->IsLeaderNode();
	//bool test_B = p_node_B->IsLeaderNode();
	//bool test = test_A + test_B;
	//std::cout<<"test\t"<<nodes<<std::endl;
	
	//std::cout<<"cell_size\t"<<cell_size<<"\t"<<"cell_threshold_division_length\t"<<cell_iter->GetThresholdDivisionLength()<<std::endl;
	
	if (cell_size >= cell_iter->GetThresholdDivisionLength() && isleader == false)
	{
		// Store parent ID for output if required
                //unsigned parent_cell_id = cell_iter->GetCellId();
		std::cout<<"inside cell_size >= "<<std::endl;
                // Create a new cell
                CellPtr p_new_cell = cell_iter->Divide_ab();

		// Add the new cell to the cell population
                mrCellPopulation.AddCell_ab(p_new_cell, *cell_iter);

		// Update counter
                num_births_this_step++;

	}
	
        
    }
	//exit(0);
    return num_births_this_step;
}





template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::DoCellBirth()
{
    	//std::cout<<"inside_docellbirth"<<std::endl;
    if (mNoBirth)
    {
        return 0;
    }

    unsigned num_births_this_step = 0;

    // Iterate over all cells, seeing if each one can be divided
    for (typename AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
	//std::cout<<"inside_docellbirth_loop"<<std::endl;
        // Check if this cell is ready to divide
        double cell_age = cell_iter->GetAge();
        if (cell_age > 0.0)
        {
		//std::cout<<"inside_docellbirth_cell_age_loop"<<std::endl;
            if (cell_iter->ReadyToDivide())
            {
		//std::cout<<"inside_docellbirth_ready_to_divide_loop"<<std::endl;
                // Check if there is room into which the cell may divide
                if (mrCellPopulation.IsRoomToDivide(*cell_iter))
                {
                    // Store parent ID for output if required
                    unsigned parent_cell_id = cell_iter->GetCellId();

                    // Create a new cell
                    CellPtr p_new_cell = cell_iter->Divide();

                    /**
                     * If required, output this location to file
                     *
                     * \todo (#2578)
                     *
                     * For consistency with the rest of the output code, consider removing the
                     * AbstractCellBasedSimulation member mOutputDivisionLocations, adding a new
                     * member mAgesAndLocationsOfDividingCells to AbstractCellPopulation, adding
                     * a new class CellDivisionLocationsWriter to the CellPopulationWriter hierarchy
                     * to output the content of mAgesAndLocationsOfDividingCells to file (remembering
                     * to clear mAgesAndLocationsOfDividingCells at each timestep), and replacing the
                     * following conditional statement with something like
                     *
                     * if (mrCellPopulation.HasWriter<CellDivisionLocationsWriter>())
                     * {
                     *     mCellDivisionLocations.push_back(new_location);
                     * }
                     */
                    if (mOutputDivisionLocations)
                    {
                        c_vector<double, SPACE_DIM> cell_location = mrCellPopulation.GetLocationOfCellCentre(*cell_iter);

                        *mpDivisionLocationFile << SimulationTime::Instance()->GetTime() << "\t";
                        for (unsigned i=0; i<SPACE_DIM; i++)
                        {
                            *mpDivisionLocationFile << cell_location[i] << "\t";
                        }
                        *mpDivisionLocationFile << "\t" << cell_age << "\t" << parent_cell_id << "\t" << cell_iter->GetCellId() << "\t" << p_new_cell->GetCellId() << "\n";
                    }

                    // Add the new cell to the cell population
                    mrCellPopulation.AddCell(p_new_cell, *cell_iter);

                    // Update counter
                    num_births_this_step++;
                }
            }
        }
    }
    return num_births_this_step;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::DoCellRemoval()
{
    unsigned num_deaths_this_step = 0;

    /*
     * This labels cells as dead or apoptosing. It does not actually remove the cells,
     * mrCellPopulation.RemoveDeadCells() needs to be called for this.
     */
    for (typename std::vector<boost::shared_ptr<AbstractCellKiller<SPACE_DIM> > >::iterator killer_iter = mCellKillers.begin();
         killer_iter != mCellKillers.end();
         ++killer_iter)
    {
        (*killer_iter)->CheckAndLabelCellsForApoptosisOrDeath();
    }

    num_deaths_this_step += mrCellPopulation.RemoveDeadCells();

    return num_deaths_this_step;
}

//added 12-10-21
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::SetNumTimeSteps(unsigned numTimeSteps)
{
    //assert(dt > 0);
    mNumTimeSteps = numTimeSteps;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::SetDt(double dt)
{
    assert(dt > 0);
    mDt = dt;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::GetDt()
{
    return mDt;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::GetNumBirths()
{
    return mNumBirths;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::GetNumDeaths()
{
    return mNumDeaths;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::SetEndTime(double endTime)
{
    assert(endTime > 0);
    mEndTime = endTime;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::SetOutputDirectory(std::string outputDirectory)
{
    mOutputDirectory = outputDirectory;
    mSimulationOutputDirectory = mOutputDirectory;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::string AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::GetOutputDirectory()
{
    return mOutputDirectory;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::SetSamplingTimestepMultiple(unsigned samplingTimestepMultiple)
{
    assert(samplingTimestepMultiple > 0);
    mSamplingTimestepMultiple = samplingTimestepMultiple;
	//std::cout<<"check"<<std::endl;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::rGetCellPopulation()
{
    return mrCellPopulation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::rGetCellPopulation() const
{
    return mrCellPopulation;
}

//added 12-10-21
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::SetUpdateAdhesion1(bool updateAdhesion1)
{
    mUpdateAdhesion1 = updateAdhesion1;
}


//added 11-30-21
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::SetThresholdBasedCellDivision(bool thresholdBasedCellDivision)
{
    mthresholdBasedCellDivision = thresholdBasedCellDivision;
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::SetUpdateCellPopulationRule(bool updateCellPopulation)
{
    mUpdateCellPopulation = updateCellPopulation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::GetUpdateCellPopulationRule()
{
    return mUpdateCellPopulation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::SetNoBirth(bool noBirth)
{
    mNoBirth = noBirth;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::AddCellKiller(boost::shared_ptr<AbstractCellKiller<SPACE_DIM> > pCellKiller)
{
    mCellKillers.push_back(pCellKiller);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::RemoveAllCellKillers()
{
    mCellKillers.clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::AddSimulationModifier(boost::shared_ptr<AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM> > pSimulationModifier)
{
    mSimulationModifiers.push_back(pSimulationModifier);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM> > >* AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::GetSimulationModifiers()
{
    return &mSimulationModifiers;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::GetNodeLocation(const unsigned& rNodeIndex)
{
    std::vector<double> location;
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        location.push_back(mrCellPopulation.GetNode(rNodeIndex)->rGetLocation()[i]);
    }
    return location;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::Solve()
{
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::EVERYTHING);
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::SETUP);

    // Set up the simulation time
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetTime();

    assert(mDt != DOUBLE_UNSET);  //Subclass constructors take care of this

    if (mEndTime == DOUBLE_UNSET)
    {
        EXCEPTION("SetEndTime has not yet been called.");
    }

    /*
     * Note that mDt is used here for "ideal time step". If this step doesn't divide the time remaining
     * then a *different* time step will be taken by the time-stepper. The real time-step (used in the
     * SimulationTime singleton) is currently not available to this class.
     *
     * \todo Should we over-write the value of mDt, or change this behaviour? (see #2159)
     */
    unsigned num_time_steps = (unsigned) ((mEndTime-current_time)/mDt+0.5);
    if (current_time > 0) // use the reset function if necessary
    {
        p_simulation_time->ResetEndTimeAndNumberOfTimeSteps(mEndTime, num_time_steps);
    }
    else
    {
        if (p_simulation_time->IsEndTimeAndNumberOfTimeStepsSetUp())
        {
            EXCEPTION("End time and number of timesteps already setup. You should not use SimulationTime::SetEndTimeAndNumberOfTimeSteps in cell-based tests.");
        }
        else
        {
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(mEndTime, num_time_steps);
        }
    }

//**************SETTING UP OUTPUTDIRECTORY AND RESULTS FOLDER FOR CELL, COLLAGEN AND SUBSTRATE************************************************************

    if (mOutputDirectory == "")
    {
        EXCEPTION("OutputDirectory not set");
    }

    double time_now = p_simulation_time->GetTime();
    std::ostringstream time_string;
    time_string << time_now;
	std::string substrate_results_directory = mOutputDirectory +"/substrate_results_from_time_" + time_string.str();    //added 7-29-20
	mSubstrateSimulationOutputDirectory = substrate_results_directory; //added 7-29-20 
	OutputFileHandler output_file_handler_ab(substrate_results_directory+"/", true);
	mpVizSetupFile_ab = output_file_handler_ab.OpenOutputFile("results_ab.vizsetup");
	InitiateVtkMetaFile_ab(output_file_handler_ab); //added 7-31-20

	OutputFileHandler output_file_handler_ab2(substrate_results_directory+"/", true); 		//added 8-24-21
	InitiateSubPosFile_ab(output_file_handler_ab2); 		//added 8-24-21
	WriteSubPosToFile_ab(output_file_handler_ab2, time_now); 		//added 8-25-21

	OutputFileHandler output_file_handler_ab3(substrate_results_directory+"/", true);
	mpSubForcesFile = output_file_handler_ab3.OpenOutputFile("substratenodeforces.dat");    //added 8-27-21

	//OutputFileHandler output_file_handler_ab5(substrate_results_directory+"/", true);    //added 3-24-22
	//mpSubPositionsFile = output_file_handler_ab5.OpenOutputFile("substratenodepositions.dat");    //added 3-24-22

	//std::cout<<"inside solve_printing substrate_result_directory\t"<<substrate_results_directory<<std::endl;
	//std::cout<<"inside solve_printing mOutputDirectory\t"<<mOutputDirectory<<std::endl;
//exit(0);	
	//mpSubPositionsFile->close();
//exit(0);
	std::string collagen_results_directory = mOutputDirectory +"/collagen_results_from_time_" + time_string.str(); 		   //added 7-26-21
	mCollagenSimulationOutputDirectory = collagen_results_directory;			 //added 7-26-21 
	OutputFileHandler output_file_handler_ab1(collagen_results_directory+"/", true);
	mpVizSetupFile_ab1 = output_file_handler_ab1.OpenOutputFile("results_ab1.vizsetup");
	InitiateVtkMetaFile_ab1(output_file_handler_ab1); 		//added 7-26-21

	//OutputFileHandler output_file_handler_ab12(collagen_results_directory+"/", true); 		//added 8-24-21
	//InitiateColVelFile_ab1(output_file_handler_ab12); 		//added 8-24-21


//exit(0);

/*
#ifdef CHASTE_VTK
    mpVtkMetaFile_ab = output_file_handler_ab.OpenOutputFile("results_ab.pvd");
    *mpVtkMetaFile_ab << "<?xml version=\"1.0\"?>\n";
    *mpVtkMetaFile_ab << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
    *mpVtkMetaFile_ab << "    <Collection>\n";
#endif //CHASTE_VTK
*/
	//std::cout<<"abcd"<<std::endl;
	WriteVtkResultsToFile_ab(substrate_results_directory); //added 7-29-20

	WriteVtkResultsToFile_ab1(collagen_results_directory); 		//added 7-26-21	
	//exit(0);
	//std::cout<<"abcd2"<<std::endl;
/*
#ifdef CHASTE_VTK
    *mpVtkMetaFile_ab << "    </Collection>\n";
    *mpVtkMetaFile_ab << "</VTKFile>\n";
    mpVtkMetaFile_ab->close();
#endif //CHASTE_VTK
*/
	CloseVtkMetaFile_ab(); //added 7-31-20
	*mpVizSetupFile_ab << std::flush;


	CloseVtkMetaFile_ab1(); 		//added 7-26-21
	*mpVizSetupFile_ab1 << std::flush;
	

    std::string results_directory = mOutputDirectory +"/results_from_time_" + time_string.str();
    mSimulationOutputDirectory = results_directory;

//*****************************************************************************************************************************************************
	
	//std::cout<<"results_directory"<<results_directory<<std::endl;
    // Set up simulation

    // Create output files for the visualizer
    OutputFileHandler output_file_handler(results_directory+"/", true);

    mrCellPopulation.OpenWritersFiles(output_file_handler);   //@vertexbasedcellpopulation.cpp line 483

    if (mOutputDivisionLocations)
    {
	//std::cout<<"mOutputDivisionLocations"<<std::endl;          
	mpDivisionLocationFile = output_file_handler.OpenOutputFile("divisions.dat");// @.hpp out_stream mpDivisionLocationFile;
    }
    if (mOutputCellVelocities)
    {
        OutputFileHandler output_file_handler2(this->mSimulationOutputDirectory+"/", false);
        mpCellVelocitiesFile = output_file_handler2.OpenOutputFile("cellvelocities.dat");
    }
	//std::cout<<"abcd3"<<std::endl;
	OutputFileHandler output_file_handler3(this->mSimulationOutputDirectory+"/", false);		//added 9-12-21
	mpCellNodeForcesFile = output_file_handler3.OpenOutputFile("cellnodeforces.dat");		//added 9-12-21

	OutputFileHandler output_file_handler4(this->mSimulationOutputDirectory+"/", false);		//added 9-12-21
	mpCellNodePositionsFile = output_file_handler4.OpenOutputFile("cellnodepositions.dat");		//added 9-12-21

	unsigned index = 0;	
	for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = 	mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
             node_iter != mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
             ++node_iter,  ++index)		//added 9-12-21
        {
		//std::cout<<"abcd4"<<std::endl;
		const c_vector<double,SPACE_DIM>& position = (node_iter)->rGetLocation();
		*mpCellNodePositionsFile << position[0] << " ";		
		//std::cout<<"He-Man"<<std::endl; 
		//Node<SPACE_DIM>* p_cell_node_ab = mrCellPopulation.rGetMesh().GetNode((node_iter)->GetIndex());
		//unsigned sub_index = p_cell_node_ab->rGetSubstrateNodeIndex();
		//std::cout<<"inside abstractcell"<<std::endl;
		//std::cout<<"cell node index is\t"<<(node_iter)->GetIndex()<<"\t"<<"nearest_sub_node_index is\t"<<sub_index<<"IsFocal\t"<<p_cell_node_ab->IsFocal()<<std::endl;  
		//std::cout<<"cell node index is\t"<<(node_iter)->GetIndex()<<std::endl;       
		//old_node_locations[&(*node_iter)] = (node_iter)->rGetLocation();
        }
	*mpCellNodePositionsFile << "\n";
	*mpCellNodePositionsFile << "\n";	//added 11-23-21 to facilitate processing in matlab, having a blank line between time points will ease in post processing
    if (PetscTools::AmMaster())
    {
        mpVizSetupFile = output_file_handler.OpenOutputFile("results.vizsetup");

        for (typename std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM> > >::iterator iter = mSimulationModifiers.begin();
             iter != mSimulationModifiers.end();
             ++iter)
        {
            if (boost::dynamic_pointer_cast<AbstractPdeModifier<SPACE_DIM> >(*iter))
            {
		//std::cout<<"*this->mpVizSetupFile"<<std::endl;
                *this->mpVizSetupFile << "PDE \n";
            }
        }
    }

    this->mrCellPopulation.SimulationSetupHook(this);//@vertexbasedcellpopulation.cpp line 832


//*************************************** SETUP SOLVE AND SETUP COLLAGEN SUBSTRATE INTERACTION ******************************************************************
    SetupSolve();

	

	if(mIsCol)
	{
		SetupCellCollagenLink(); 		//added 8-4-21
		SetupCollagenSubstrateLink();
	}

	//exit(0);

//***************************************************************************************************************************************************************

    // Call SetupSolve() on each modifier
    for (typename std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM> > >::iterator iter = mSimulationModifiers.begin();
         iter != mSimulationModifiers.end();
         ++iter)
    {
        (*iter)->SetupSolve(this->mrCellPopulation,this->mSimulationOutputDirectory);
    }

    /*
     * Age the cells to the correct time. Note that cells are created with
     * negative birth times so that some are initially almost ready to divide.
     */
    LOG(1, "Setting up cells...");
    for (typename AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
        /*
         * We don't use the result; this call is just to force the cells to age
         * to the current time running their cell-cycle models to get there.
         */
        cell_iter->ReadyToDivide();
	//std::cout<<"inside ready-to-divide"<<std::endl;
    }
    LOG(1, "\tdone\n");
	//exit(0);

    // Write initial conditions to file for the visualizer
    WriteVisualizerSetupFile();  //@offlatticesimulation_ab.cpp line 335
	//std::cout<<"PetscTools::AmMaster1()"<<"\t"<<PetscTools::AmMaster()<<std::endl;
    if (PetscTools::AmMaster())
    {
        *mpVizSetupFile << std::flush;
    }
	//std::cout<<"PetscTools::AmMaster2()"<<"\t"<<PetscTools::AmMaster()<<std::endl;
    mrCellPopulation.WriteResultsToFiles(results_directory+"/");   //abstractcellpopulation.cpp line 589

    OutputSimulationSetup();
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::SETUP);




//****************************************************************************************************************************************************************
//*************************************************   ENTER MAIN TIME LOOP  *************************************************************************************
//****************************************************************************************************************************************************************

unsigned abc = 1;
    while (!( p_simulation_time->IsFinished() || StoppingEventHasOccurred() ) )
    {
        LOG(1, "--TIME = " << p_simulation_time->GetTime() << "\n");
	//std::cout<<"tipsyk1"<<std::endl;

	


        // This function calls DoCellRemoval(), DoCellBirth() and CellPopulation::Update()
        //UpdateCellPopulation();
	//UpdateCellPopulation_ab();

	if (mthresholdBasedCellDivision)
	{
		UpdateCellPopulation_ab();
		SetupCellCollagenLink();
	}
	
	
	//std::cout<<"tipsyk2"<<std::endl;
	//std::cout<<"after_UpdateCellPopulation"<<std::endl;

/*
	//added 8-7-20
	for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = 	mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
             node_iter != mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
        {
		//std::cout<<"He-Man"<<std::endl; 
		//Node<SPACE_DIM>* p_cell_node_ab = mrCellPopulation.rGetMesh().GetNode((node_iter)->GetIndex());
		//unsigned sub_index = p_cell_node_ab->rGetSubstrateNodeIndex();
		//std::cout<<"inside abstractcell"<<std::endl;
		//std::cout<<"cell node index is\t"<<(node_iter)->GetIndex()<<"\t"<<"nearest_sub_node_index is\t"<<sub_index<<"IsFocal\t"<<p_cell_node_ab->IsFocal()<<std::endl;  
		std::cout<<"cell node index is\t"<<(node_iter)->GetIndex()<<std::endl;       
		//old_node_locations[&(*node_iter)] = (node_iter)->rGetLocation();
        }
//exit(0);
*/




	//if(mIsCol)		//commented 12-10-21
	if(mUpdateAdhesion1)		//added 12-10-21
	{
		//std::cout<<"tipsyk2.1"<<std::endl;	
		//SetupCellCollagenLink(); 		//added 9-24-21	
		UpdateAdhesions1();     //added 7-15-21
		//std::cout<<"tipsyk2.5"<<std::endl;
		//exit(0);
	}
	/*		//commented 12-10-21
	else
	{
		UpdateAdhesions();     //added 7-14-20
	}*/
	//std::cout<<"tipsyk33"<<std::endl;
//exit(0);
        // Store whether we are sampling results at the current timestep
        SimulationTime* p_time = SimulationTime::Instance();
        bool at_sampling_timestep = (p_time->GetTimeStepsElapsed()%this->mSamplingTimestepMultiple == 0);

        /*
         * If required, store the current locations of cell centres. Note that we need to
         * use a std::map between cells and locations, rather than (say) a std::vector with
         * location indices corresponding to cells, since once we call UpdateCellLocations()
         * the location index of each cell may change. This is especially true in the case
         * of a CaBasedCellPopulation.
         */
        std::map<CellPtr, c_vector<double, SPACE_DIM> > old_cell_locations;
	std::map<CellPtr, c_vector<double, SPACE_DIM> > old_collagen_locations;
	std::map<CellPtr, c_vector<double, SPACE_DIM> > old_substrate_locations;
	

	if (mIsCol)
	{

		if (mOutputCellVelocities && at_sampling_timestep)
        	{
            		//std::cout<<"inside moutputcellvelocities"<<std::endl;
	    		for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = mrCellPopulation.Begin();
                 cell_iter != mrCellPopulation.End();
                 ++cell_iter)
            		{
                		old_cell_locations[*cell_iter] = mrCellPopulation.GetLocationOfCellCentre_ab(*cell_iter);
				//std::cout<<"old_cell_locations[*cell_iter]\t"<<old_cell_locations[*cell_iter]<<std::endl;
           		}

       		}
	}

	
	else
	{
        	if (mOutputCellVelocities && at_sampling_timestep)
        	{
            		for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = mrCellPopulation.Begin();
                 cell_iter != mrCellPopulation.End();
                 ++cell_iter)
            		{
                		old_cell_locations[*cell_iter] = mrCellPopulation.GetLocationOfCellCentre(*cell_iter);
            		}
        	}
	}

//exit(0);

	//std::cout<<"He-Man2"<<std::endl;


        // Update cell locations and topology
	//std::cout<<"tipsyk4"<<std::endl;
        UpdateCellLocationsAndTopology();
	//exit(0);
	
	///std::cout<<"Inside Solve() in AbstractCellBasedSimulation tipsyk4"<<std::endl;
        UpdateCollagenTopology();


	//std::cout<<"tipsyk5"<<std::endl;
	//exit(0);
	UpdateSubstrateTopology();   //added 8-2-20
	
	//std::cout<<"tipsyk5.3"<<std::endl;	
	
	//exit(0);
	//std::cout<<"after_UpdateSubstrateTopology()"<<std::endl;



//exit(0);
	//UpdateCellLocationsAndTopology_ab();    //added 6-17-20

        // Now write cell velocities to file if required
	//std::cout<<mOutputCellVelocities<<std::endl;
        if (mOutputCellVelocities && at_sampling_timestep)
        {
            time_now = p_time->GetTime() + mDt;		//added 8-26-21
		//std::cout<<"mDt\t"<<mDt<<std::endl;
		//exit(0);
	    WriteSubPosToFile_ab(output_file_handler_ab2, time_now);		//added 8-26-21
	    // Offset as doing this before we increase time by mDt
            *mpCellVelocitiesFile << p_time->GetTime() + mDt<< "\t";
		//std::cout<<"mpCellVelocitiesFile"<<std::endl;
            for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = mrCellPopulation.Begin();
                 cell_iter != mrCellPopulation.End();
                 ++cell_iter)
            {
                unsigned index = mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                const c_vector<double,SPACE_DIM>& position = mrCellPopulation.GetLocationOfCellCentre_ab(*cell_iter);

                c_vector<double, SPACE_DIM> velocity; // Two lines for profile build
                velocity = (position - old_cell_locations[*cell_iter])/mDt;

                *mpCellVelocitiesFile << index  << " ";
                for (unsigned i=0; i<SPACE_DIM; i++)
                {
                    *mpCellVelocitiesFile << position[i] << " ";

			//std::cout<<"inside mpcellvelocitiesfile"<<std::endl;

                }

                for (unsigned i=0; i<SPACE_DIM; i++)
                {
                    *mpCellVelocitiesFile << velocity[i] << " ";
                }
            }
            *mpCellVelocitiesFile << "\n";
        }

	//std::cout<<"tipsyk5.4"<<std::endl;

	for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = 	mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
             node_iter != mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
             ++node_iter,  ++index)		//added 9-13-21
        {
		
		const c_vector<double,SPACE_DIM>& position = (node_iter)->rGetLocation();
		*mpCellNodePositionsFile << position[0] << " ";		
		//std::cout<<"He-Man"<<std::endl; 
		//Node<SPACE_DIM>* p_cell_node_ab = mrCellPopulation.rGetMesh().GetNode((node_iter)->GetIndex());
		//unsigned sub_index = p_cell_node_ab->rGetSubstrateNodeIndex();
		//std::cout<<"inside abstractcell"<<std::endl;
		//std::cout<<"cell node index is\t"<<(node_iter)->GetIndex()<<"\t"<<"nearest_sub_node_index is\t"<<sub_index<<"IsFocal\t"<<p_cell_node_ab->IsFocal()<<std::endl;  
		//std::cout<<"cell node index is\t"<<(node_iter)->GetIndex()<<std::endl;       
		//old_node_locations[&(*node_iter)] = (node_iter)->rGetLocation();
        }
	*mpCellNodePositionsFile << "\n";
	*mpCellNodePositionsFile << "\n";	//added 11-23-21 to facilitate processing in matlab, having a blank line between time points will ease in post processing

	//std::cout<<"outside mpcellvelocitiesfile"<<std::endl;

	//exit(0);
        // Update the assignment of cells to processes.
        mrCellPopulation.UpdateCellProcessLocation();
	//std::cout<<"after_UpdateCellProcessLocation"<<std::endl;

	
        // Increment simulation time here, so results files look sensible
        p_simulation_time->IncrementTimeOneStep();

        // Call UpdateAtEndOfTimeStep() on each modifier
        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::UPDATESIMULATION);
        for (typename std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM> > >::iterator iter = mSimulationModifiers.begin();
             iter != mSimulationModifiers.end();
             ++iter)
        {
		//std::cout<<"b4_UpdateAtEndOfTimeStep"<<std::endl;
            (*iter)->UpdateAtEndOfTimeStep(this->mrCellPopulation);
		//std::cout<<"after_UpdateAtEndOfTimeStep"<<std::endl;
        }
	//std::cout<<"b4_EndEvent"<<std::endl;
        CellBasedEventHandler::EndEvent(CellBasedEventHandler::UPDATESIMULATION);
	//std::cout<<"after_EndEvent"<<std::endl;
        // Output current results to file
        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::OUTPUT);
	//std::cout<<"after_BeginEvent"<<std::endl;
	//std::cout<<p_simulation_time->GetTimeStepsElapsed()%mSamplingTimestepMultiple<<std::endl;
	//std::cout<<"this->mSamplingTimestepMultiple\t"<<this->mSamplingTimestepMultiple<<"\t"<<"mSamplingTimestepMultiple\t"<<mSamplingTimestepMultiple<<std::endl;
        if (p_simulation_time->GetTimeStepsElapsed()%mSamplingTimestepMultiple == 0)// should be at_sampling_timestep !
        {
		//std::cout<<"b4_WriteResultsToFiles"<<std::endl;
            mrCellPopulation.WriteResultsToFiles(results_directory+"/"); //abstractcellpopulation.cpp line 589
		//std::cout<<"after_WriteResultsToFiles"<<std::endl;
		
		WriteVtkResultsToFile_ab(substrate_results_directory);   //added 7-31-20

		WriteVtkResultsToFile_ab1(collagen_results_directory); 		  //added 7-26-21
	
		//exit(0);

		//std::cout<<"after_WriteVtkResultsToFile_ab"<<std::endl;	
            // Call UpdateAtEndOfOutputTimeStep() on each modifier
            for (typename std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM> > >::iterator iter = mSimulationModifiers.begin();
                 iter != mSimulationModifiers.end();
                 ++iter)
            {
                (*iter)->UpdateAtEndOfOutputTimeStep(this->mrCellPopulation);
            }
        }
	//std::cout<<"b4_EndEvent"<<std::endl;
        CellBasedEventHandler::EndEvent(CellBasedEventHandler::OUTPUT);
	//std::cout<<"after_EndEvent"<<std::endl;

	//std::cout<<"abc\t"<<abc<<std::endl;

///*
	if (abc == mNumTimeSteps)
	{
		exit(0);
	}
	abc++;
//*/


    }
	//std::cout<<"b4_Endtime"<<std::endl;
    LOG(1, "--END TIME = " << p_simulation_time->GetTime() << "\n");
	//std::cout<<"after_Endtime"<<std::endl;
    /*
     * Carry out a final update so that cell population is coherent with new cell positions.
     * Note that cell birth and death still need to be checked because they may be spatially
     * dependent.
     */
	//std::cout<<"b4_UpdateCellPopulation"<<std::endl;
    //UpdateCellPopulation();
	//std::cout<<"after_UpdateCellPopulation"<<std::endl;




    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::UPDATESIMULATION);
    // Call UpdateAtEndOfSolve(), on each modifier
    for (typename std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM> > >::iterator iter = mSimulationModifiers.begin();
         iter != mSimulationModifiers.end();
         ++iter)
    {
        (*iter)->UpdateAtEndOfSolve(this->mrCellPopulation);
    }
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::UPDATESIMULATION);

    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::OUTPUT);

    mrCellPopulation.CloseWritersFiles();

    if (mOutputDivisionLocations)
    {
        mpDivisionLocationFile->close();
    }
    if (mOutputCellVelocities)
    {
        mpCellVelocitiesFile->close();
    }

	mpSubPositionsFile->close(); 	//added 8-26-21

	mpSubForcesFile->close(); 	//added 8-27-21

	mpCellNodeForcesFile->close(); 	//added 9-13-21

	mpCellNodePositionsFile->close();	//added 9-13-21

    if (PetscTools::AmMaster())
    {
        *mpVizSetupFile << "Complete\n";
        mpVizSetupFile->close();
    }
	CloseVtkMetaFile_ab(); //added 7-31-20

	CloseVtkMetaFile_ab1(); 		//added 7-26-21

	//std::cout<<"CloseVtkMetaFile_ab()"<<std::endl;
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::OUTPUT);
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::EVERYTHING);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::StoppingEventHasOccurred()
{
    return false;
}
/*
//added 7-14-20
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::UpdateAdhesions()
{
	std::cout<<"nemesis"<<std::endl;
}
*/

//added 9-21-21
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::UpdateCellPopulation_ab()
{
	//std::cout<<"inside updatecellpopulation_ab()"<<std::endl;
    /*// Remove dead cells
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::DEATH);
    unsigned deaths_this_step = DoCellRemoval();
    mNumDeaths += deaths_this_step;
    LOG(1, "\tNum deaths = " << mNumDeaths << "\n");
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::DEATH);*/

    // Divide cells
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::BIRTH);
    unsigned births_this_step = DoCellBirth_ab();
	//std::cout<<"inside_updatecellpopulation_ab_printing_births_this_step\t"<<births_this_step<<std::endl;
    mNumBirths += births_this_step;
    LOG(1, "\tNum births = " << mNumBirths << "\n");
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::BIRTH);

    // This allows NodeBasedCellPopulation::Update() to do the minimum amount of work
    //bool births_or_death_occurred = ((births_this_step>0) || (deaths_this_step>0));

	bool births_or_death_occurred = (births_this_step>0);	//added 9-21-21
	//std::cout<<"inside_updatecellpopulation_before_update_loop"<<std::endl;
    // Update topology of cell population
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::UPDATECELLPOPULATION);
    if (mUpdateCellPopulation)
    {
        LOG(1, "\tUpdating cell population...");
        mrCellPopulation.Update_ab(births_or_death_occurred);
        LOG(1, "\tdone.\n");
    }
    else if (births_or_death_occurred)
    {
        EXCEPTION("CellPopulation has had births or deaths but mUpdateCellPopulation is set to false, please set it to true.");
    }
	//std::cout<<"inside_updatecellpopulation_after_update_loop"<<std::endl;
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::UPDATECELLPOPULATION);
}





template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::UpdateCellPopulation()
{
	std::cout<<"inside updatecellpopulation()"<<std::endl;
    // Remove dead cells
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::DEATH);
    unsigned deaths_this_step = DoCellRemoval();
    mNumDeaths += deaths_this_step;
    LOG(1, "\tNum deaths = " << mNumDeaths << "\n");
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::DEATH);

    // Divide cells
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::BIRTH);
    unsigned births_this_step = DoCellBirth();
    mNumBirths += births_this_step;
    LOG(1, "\tNum births = " << mNumBirths << "\n");
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::BIRTH);

    // This allows NodeBasedCellPopulation::Update() to do the minimum amount of work
    bool births_or_death_occurred = ((births_this_step>0) || (deaths_this_step>0));

    // Update topology of cell population
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::UPDATECELLPOPULATION);
    if (mUpdateCellPopulation)
    {
        LOG(1, "\tUpdating cell population...");
        mrCellPopulation.Update(births_or_death_occurred);
        LOG(1, "\tdone.\n");
    }
    else if (births_or_death_occurred)
    {
        EXCEPTION("CellPopulation has had births or deaths but mUpdateCellPopulation is set to false, please set it to true.");
    }
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::UPDATECELLPOPULATION);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::GetOutputDivisionLocations()
{
    return mOutputDivisionLocations;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::SetOutputDivisionLocations(bool outputDivisionLocations)
{
    mOutputDivisionLocations = outputDivisionLocations;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::GetOutputCellVelocities()
{
    return mOutputCellVelocities;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::SetOutputCellVelocities(bool outputCellVelocities)
{
    mOutputCellVelocities = outputCellVelocities;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::OutputSimulationSetup()
{
    OutputFileHandler output_file_handler(this->mSimulationOutputDirectory + "/", false);

    // Output machine information
    ExecutableSupport::SetOutputDirectory(output_file_handler.GetOutputDirectoryFullPath());
    ExecutableSupport::WriteMachineInfoFile("system_info");

    if (PetscTools::AmMaster())
    {
        // Output Chaste provenance information
        out_stream build_info_file = output_file_handler.OpenOutputFile("build.info");
        std::string build_info;
        ExecutableSupport::GetBuildInfo(build_info);
        *build_info_file << build_info;
        build_info_file->close();

        // Output simulation parameter and setup details
        out_stream parameter_file = output_file_handler.OpenOutputFile("results.parameters");

        // Output simulation details
        std::string simulation_type = GetIdentifier();

        *parameter_file << "<Chaste>\n";
        *parameter_file << "\n\t<" << simulation_type << ">\n";
        OutputSimulationParameters(parameter_file);
        *parameter_file << "\t</" << simulation_type << ">\n";
        *parameter_file << "\n";

        // Output cell population details (includes cell-cycle model details)
        mrCellPopulation.OutputCellPopulationInfo(parameter_file);

        // Loop over cell killers
        *parameter_file << "\n\t<CellKillers>\n";
        for (typename std::vector<boost::shared_ptr<AbstractCellKiller<SPACE_DIM> > >::iterator iter = mCellKillers.begin();
             iter != mCellKillers.end();
             ++iter)
        {
            // Output cell killer details
            (*iter)->OutputCellKillerInfo(parameter_file);
        }
        *parameter_file << "\t</CellKillers>\n";

        // Iterate over simulationmodifiers
        *parameter_file << "\n\t<SimulationModifiers>\n";
        for (typename std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM> > >::iterator iter = mSimulationModifiers.begin();
             iter != mSimulationModifiers.end();
             ++iter)
        {
            // Output simulation modifier details
            (*iter)->OutputSimulationModifierInfo(parameter_file);
        }
        *parameter_file << "\t</SimulationModifiers>\n";

        // This is used to output information about subclasses
        OutputAdditionalSimulationSetup(parameter_file);

        *parameter_file << "\n</Chaste>\n";
        parameter_file->close();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::OutputSimulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<Dt>" << mDt << "</Dt>\n";
    *rParamsFile << "\t\t<EndTime>" << mEndTime << "</EndTime>\n";
    *rParamsFile << "\t\t<SamplingTimestepMultiple>" << mSamplingTimestepMultiple << "</SamplingTimestepMultiple>\n";
    *rParamsFile << "\t\t<OutputDivisionLocations>" << mOutputDivisionLocations << "</OutputDivisionLocations>\n";
    *rParamsFile << "\t\t<OutputCellVelocities>" << mOutputCellVelocities << "</OutputCellVelocities>\n";
}

// Explicit instantiation
template class AbstractCellBasedSimulation<1,1>;
template class AbstractCellBasedSimulation<1,2>;
template class AbstractCellBasedSimulation<2,2>;
template class AbstractCellBasedSimulation<1,3>;
template class AbstractCellBasedSimulation<2,3>;
template class AbstractCellBasedSimulation<3,3>;
