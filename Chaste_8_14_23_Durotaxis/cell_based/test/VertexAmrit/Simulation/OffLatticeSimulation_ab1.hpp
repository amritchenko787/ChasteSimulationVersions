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

#ifndef OFFLATTICESIMULATION_AB1_HPP_
#define OFFLATTICESIMULATION_AB1_HPP_

#include "MutableVertexMesh.hpp"
#include "AbstractCellBasedSimulation.hpp"
#include "AbstractForce.hpp"
#include "AbstractCellPopulationBoundaryCondition.hpp"
#include "AbstractNumericalMethod.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>
#include "UblasVectorInclude.hpp" //added 7-14-20
#include "UblasMatrixInclude.hpp" //added 7-14-20
#include <boost/utility.hpp> //added 7-14-20

/**
 * Run an off-lattice 2D or 3D cell-based simulation using an off-lattice
 * cell population.
 *
 * In cell-centre-based cell populations, each cell is represented by a
 * single node (corresponding to its centre), and connectivity is defined
 * either by a Delaunay triangulation or a radius of influence. In vertex-
 * based cell populations, each cell is represented by a polytope
 * (corresponding to its membrane) with a variable number of vertices.
 * Alternative cell populations may be defined by the user.
 *
 * The OffLatticeSimulation_ab1 is constructed with a CellPopulation, which
 * updates the correspondence between each Cell and its spatial representation
 * and handles cell division (governed by the CellCycleModel associated
 * with each cell). Once constructed, one or more Force laws may be passed
 * to the OffLatticeSimulation_ab1 object, to define the mechanical properties
 * of the CellPopulation. Similarly, one or more CellKillers may be passed
 * to the OffLatticeSimulation_ab1 object to specify conditions in which Cells
 * may die, and one or more CellPopulationBoundaryConditions to specify
 * regions in space beyond which Cells may not move.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class OffLatticeSimulation_ab1 : public AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    friend class TestOffLatticeSimulation_ab1;
    friend class TestOffLatticeSimulation_ab1WithNodeBasedCellPopulation;

    /**
     * Save or restore the simulation.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM> >(*this);
        archive & mForceCollection;
	archive & mForceCollection_ab; //added 8-4-20
	archive & mForceCollection_ab1; //added 7-14-21
        archive & mBoundaryConditions;
        archive & mpNumericalMethod;
	//archive & mpNumericalMethod_ab;//added 8-4-20
    }

protected:

    /** The mechanics used to determine the new location of the cells, a list of the forces. */
    std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > > mForceCollection;
	std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > > mForceCollection_ab; //added 8-3-20   force collection for substrate

	std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > > mForceCollection_ab1; //added 7-14-21  force collection for collagen

    /** List of boundary conditions. */
    std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> > > mBoundaryConditions;

    /** The numerical method to use in this simulation. Defaults to the explicit forward Euler method. */
    boost::shared_ptr<AbstractNumericalMethod<ELEMENT_DIM, SPACE_DIM> > mpNumericalMethod;
	boost::shared_ptr<AbstractNumericalMethod<ELEMENT_DIM, SPACE_DIM> > mpNumericalMethod_ab; //added 8-4-20
	boost::shared_ptr<AbstractNumericalMethod<ELEMENT_DIM, SPACE_DIM> > mpNumericalMethod_ab1; //added 7-14-21

    /**
     * Overridden UpdateCellLocationsAndTopology() method.
     *
     * Calculate forces and update node positions.
     */
    virtual void UpdateCellLocationsAndTopology();

	virtual void UpdateAdhesions1();     //added 7-15-21

	virtual void UpdateAdhesions();     //added 7-14-20


	virtual void UpdateCollagenTopology();     //added 7-16-21
	

	virtual void UpdateSubstrateTopology();     //added 8-2-20


 //   virtual void UpdateCellLocationsAndTopology_ab();    //added 6-17-20
    /**
     * Sends nodes back to the positions given in the input map. Used after a failed step
     * when adaptivity is turned on.
     *
     * @param oldNodeLoctions A map linking nodes to their old positions.
     */
    void RevertToOldLocations(std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > oldNodeLoctions);

    /**
     * Applies any boundary conditions.
     *
     * @param oldNodeLoctions Mapping between node indices and old node locations
     */
    void ApplyBoundaries(std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > oldNodeLoctions);

    /**
     * Overridden SetupSolve() method to clear the forces applied to the nodes.
     */
    virtual void SetupSolve();

    /**
     * Overridden WriteVisualizerSetupFile() method.
     */
    virtual void WriteVisualizerSetupFile();

    /** Reference to the substrate mesh. */      //added 6-11-20
    MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>& mrMesh_ab;

    /** Reference to the collagen mesh. */      //added 7-8-21
    MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>& mrMesh_ab1;

	virtual void InitiateVtkMetaFile_ab(OutputFileHandler& routput_file_handler_ab);    //added 7-31-20
	virtual void WriteVtkResultsToFile_ab(const std::string& rDirectory);    //added 7-29-20
 	virtual void CloseVtkMetaFile_ab();    //added 7-31-20


	virtual void InitiateVtkMetaFile_ab1(OutputFileHandler& routput_file_handler_ab1);	    //added 7-26-21
	virtual void WriteVtkResultsToFile_ab1(const std::string& rDirectory1);		    //added 7-26-21
 	virtual void CloseVtkMetaFile_ab1(); 		   //added 7-26-21
	virtual void InitiateColVelFile_ab1(OutputFileHandler& routput_file_handler_ab12);		//added 8-24-21
	virtual void InitiateSubPosFile_ab(OutputFileHandler& routput_file_handler_ab2);		//added 8-24-21
	virtual void WriteSubPosToFile_ab(OutputFileHandler& routput_file_handler_ab2, double timeNow);		//added 8-25-21

	virtual void SetupCollagenSubstrateLink();		//added 7-27-21

	virtual void SetupCellCollagenLink();		//added 8-3-21

	

public:

    /**
     * Constructor.
     *
     * @param rMesh reference to a substrate mesh (added 6//8/20)
     * @param rCellPopulation Reference to a cell population object
     * @param deleteCellPopulationInDestructor Whether to delete the cell population on destruction to
     *     free up memory (defaults to false)
     * @param initialiseCells Whether to initialise cells (defaults to true, set to false when loading
     *     from an archive)
     */
    OffLatticeSimulation_ab1(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
			 MutableVertexMesh<ELEMENT_DIM,SPACE_DIM>& rMesh_ab,
			 MutableVertexMesh<ELEMENT_DIM,SPACE_DIM>& rMesh_ab1,
                         bool deleteCellPopulationInDestructor=false,
                         bool initialiseCells=true);   // modified from the original contructor on 7-8-21, included the rMesh_ab1, supposedly the collagen mesh

    /**
     * Add a force to be used in this simulation (use this to set the mechanics system).
     *
     * @param pForce pointer to a force law
     */
    void AddForce(boost::shared_ptr<AbstractForce<ELEMENT_DIM,SPACE_DIM> > pForce);
	void AddForce_ab(boost::shared_ptr<AbstractForce<ELEMENT_DIM,SPACE_DIM> > pForce_ab);//added 7-14-20

	void AddForce_ab1(boost::shared_ptr<AbstractForce<ELEMENT_DIM,SPACE_DIM> > pForce_ab1);//added 7-14-21


	c_vector<double, SPACE_DIM> GetVectorFromAtoB_ab(const c_vector<double, SPACE_DIM>& rLocationA,
                                                          const c_vector<double, SPACE_DIM>& rLocationB);//added 7-14-20


	c_vector<double, SPACE_DIM> GetVectorFromAtoB_ab1(const c_vector<double, SPACE_DIM>& rLocationA,
                                                          const c_vector<double, SPACE_DIM>& rLocationB);	//added 7-28-21

	
	unsigned GetNearestSubstrateNode2CellNode(const c_vector<double, SPACE_DIM>& rLocationA, unsigned& rcell_node_index, bool isCol, bool isLeader); //added 7-14-20, modified 7-16-21 by adding bool isCol, modified 12-1-21 by adding bool isLeader

	
	unsigned GetNearestSubstrateNode2ColNode(const c_vector<double, SPACE_DIM>& rLocationA, unsigned& rcell_node_index);	 //added 7-27-21






    /**
     * Remove all the forces.
     */
    void RemoveAllForces();

    /**
     * Add a cell population boundary condition to be used in this simulation.
     *
     * @param pBoundaryCondition pointer to a boundary condition
     */
    void AddCellPopulationBoundaryCondition(boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> >  pBoundaryCondition);

    /**
     * Method to remove all the cell population boundary conditions
     */
    void RemoveAllCellPopulationBoundaryConditions();

    /**
     * Set the numerical method to be used in this simulation (use this to solve the mechanics system).
     *
     * @param pNumericalMethod pointer to a numerical method object
     */
    void SetNumericalMethod(boost::shared_ptr<AbstractNumericalMethod<ELEMENT_DIM, SPACE_DIM> > pNumericalMethod);

    /**
     * @return the current numerical method.
     */
    const boost::shared_ptr<AbstractNumericalMethod<ELEMENT_DIM, SPACE_DIM> > GetNumericalMethod() const;

    /**
     * Overridden OutputAdditionalSimulationSetup() method.
     *
     * Output any force, boundary condition or numerical method information.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputAdditionalSimulationSetup(out_stream& rParamsFile);

    /**
     * Overridden OutputSimulationParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputSimulationParameters(out_stream& rParamsFile);

    /**
     * Directly access the forces attached to this simulation, to allow their manipulation after archiving.
     *
     * @return mForceCollection the vector of pointers to forces attached to this simulation
     */
    const std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >& rGetForceCollection() const;
    
    MutableVertexMesh<ELEMENT_DIM,SPACE_DIM>& rGetMesh_ab();   //added 6-11-20

    const MutableVertexMesh<ELEMENT_DIM,SPACE_DIM>& rGetMesh_ab() const;  ////added 6-11-20
};

// Serialization for Boost >= 1.36
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(OffLatticeSimulation_ab1)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct an OffLatticeSimulation_ab1.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const OffLatticeSimulation_ab1<ELEMENT_DIM,SPACE_DIM> * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
    
    const MutableVertexMesh<ELEMENT_DIM,SPACE_DIM>* p_mesh_ab = &(t->rGetMesh_ab());    //added 6-11-20
    ar & p_mesh_ab;

    const MutableVertexMesh<ELEMENT_DIM,SPACE_DIM>* p_mesh_ab1 = &(t->rGetMesh_ab());    //added 7-8-21
    ar & p_mesh_ab1;	
}

/**
 * De-serialize constructor parameters and initialise an OffLatticeSimulation_ab1.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, OffLatticeSimulation_ab1<ELEMENT_DIM,SPACE_DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_cell_population;
    ar >> p_cell_population;

    MutableVertexMesh<ELEMENT_DIM,SPACE_DIM>* p_mesh_ab;         //added 6-11-20
    ar >> p_mesh_ab;

    MutableVertexMesh<ELEMENT_DIM,SPACE_DIM>* p_mesh_ab1;         //added 7-8-21
    ar >> p_mesh_ab1;	


    // Invoke inplace constructor to initialise instance, middle two variables set extra
    // member variables to be deleted as they are loaded from archive and to not initialise cells.
    ::new(t)OffLatticeSimulation_ab1<ELEMENT_DIM,SPACE_DIM>(*p_cell_population, *p_mesh_ab, *p_mesh_ab1, true, false);         //added 7-8-21
}
}
} // namespace

#endif /*OFFLATTICESIMULATION_AB1_HPP_*/
