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

#include "StraightVertexMeshGenerator_ab2.hpp"

StraightVertexMeshGenerator_ab2::StraightVertexMeshGenerator_ab2(unsigned numElementsAcross,
                                                           unsigned numElementsUp,
							   double ycoordinate,     //added 7-12-21
							   bool isSubstrate,     //added 9-16-21	
							   bool isCollagen,     //added 7-12-21, modified 8-12-21	
							   bool isLongCollagen,		//added 8-12-21	
							   double numNodesRemoved,	//added 8-13-21	
						           unsigned shortFiberConnectivity,	//added  12-25-21
                                                           bool isFlatBottom,
                                                           double cellRearrangementThreshold,
                                                           double t2Threshold,
                                                           double elementArea)
{
    assert(numElementsAcross > 0);
	assert(numElementsAcross % 4 == 0);    //added 7-1-21

    
    std::vector<Node<2>*> nodes;
    std::vector<VertexElement<2,2>*>  elements;
	std::vector<Node<2>*> element_nodes1;    //added 8-9-20
	std::vector<Node<2>*> element_nodes2;    //added 8-9-20


    unsigned node_index = 0;
    //unsigned node_indices[18];    //mod 7-1-21   
	unsigned node_indices[2];    //added 7-1-21
    unsigned element_index;
    //std::cout<<"damnation_2"<<std::endl;
    // Create the nodes, row by row, from the bottom up
	
	double kk = 0;
	double numelementsacross = numElementsAcross;

//***************************************CREATING NODES***************************************************************************************
//*********************************************************************************************************************************************

    // On the first row we have numElementsAcross nodes, all of which are boundary nodes
    for (unsigned i=0; i<numElementsAcross + 1; i++)
    {
        //if (isLongCollagen == true)
	if (isCollagen == true)		//modified 8-16-21, implementing condition for both long and short collagen
	{
		if ((kk - (numelementsacross/4)) < -(numNodesRemoved/4) || (kk - (numelementsacross/4)) > (numNodesRemoved/4))		//(kk - (numelementsacross/4)) is x-coord of node
		{
			//std::cout<<"kk is\t"<<kk<<std::endl;
			if ((kk - (numelementsacross/4)) == -((numNodesRemoved + 1)/4) || (kk - (numelementsacross/4)) == (numelementsacross)/4)	//condition for node left of the node removed and the right most node 
			{
								
				unsigned pp = 0;				
				Node<2>* p_node = new Node<2>(node_index, true, kk + (pp*0.5/3) - (numelementsacross/4), ycoordinate);
				p_node->HasFocal(false);
				p_node->HasFocal1(false);	//added 7-27-21
				p_node->HasLeaderNode(false);	//added 9-16-21
				//std::cout<<p_node->IsLeaderNode()<<std::endl;;			//added 9-16-21
        			nodes.push_back(p_node);
				//std::cout<<"collagen node index is\t"<<node_index<<"\tx-coord of node is\t"<<kk + (pp*0.5/3) - (numelementsacross/4)<<"\t"<<"y-coord of node is\t"<<ycoordinate<<std::endl;
				node_index++;
			}
			else
			{	
				for (unsigned pp=0; pp<3; pp++)
				{				
					Node<2>* p_node = new Node<2>(node_index, true, kk + (pp*0.5/3) - (numelementsacross/4), ycoordinate);
					p_node->HasFocal(false);
					p_node->HasFocal1(false);	//added 7-27-21
					p_node->HasLeaderNode(false);	//added 9-16-21
					//std::cout<<p_node->IsLeaderNode()<<std::endl;;			//added 9-16-21
        				nodes.push_back(p_node);
					//std::cout<<"collagen node index is\t"<<node_index<<"\tx-coord of node is\t"<<kk + (pp*0.5/3) - (numelementsacross/4)<<"\t"<<"y-coord of node is\t"<<ycoordinate<<std::endl;
					node_index++;
				}
			}
		}
		kk = kk + 0.5;
		//exit(0);
		
	}
	else if (isSubstrate == true)	//modified 9-16-21 substrate mesh node accounting
	{
		if ((kk - (numelementsacross/4)) == (numelementsacross)/4)
		{
			std::cout<<"kk - (numelementsacross/4)\t"<< kk - (numelementsacross/4)<<std::endl;
			unsigned pp=0;
			Node<2>* p_node = new Node<2>(node_index, true, kk + (pp*0.5/3) - (numelementsacross/4), ycoordinate);
			p_node->HasFocal(false);
			p_node->HasFocal1(false);	//added 7-27-21
			p_node->HasLeaderNode(false);	//added 9-16-21
			//std::cout<<p_node->IsLeaderNode()<<std::endl;;			//added 9-16-21
        		nodes.push_back(p_node);
			std::cout<<"substrate node index is\t"<<node_index<<"\tx-coord of node is\t"<<kk + (pp*0.5/3) - (numelementsacross/4)<<"\t"<<"y-coord of node is\t"<<ycoordinate<<std::endl;
			node_index++;
		}
		else
		{		
			for (unsigned pp=0; pp<3; pp++)
			{				
				Node<2>* p_node = new Node<2>(node_index, true, kk + (pp*0.5/3) - (numelementsacross/4), ycoordinate);
				p_node->HasFocal(false);
				p_node->HasFocal1(false);	//added 7-27-21
				p_node->HasLeaderNode(false);	//added 9-16-21
				//std::cout<<p_node->IsLeaderNode()<<std::endl;;			//added 9-16-21
        			nodes.push_back(p_node);
				std::cout<<"substrate node index is\t"<<node_index<<"\tx-coord of node is\t"<<kk + (pp*0.5/3) - (numelementsacross/4)<<"\t"<<"y-coord of node is\t"<<ycoordinate<<std::endl;
				node_index++;
			}
		}
		//Node<2>* p_node = new Node<2>(node_index, true, kk - (numelementsacross/4), ycoordinate);
		//p_node->HasFocal(false);
		//p_node->HasFocal1(false);	//added 7-27-21
		//p_node->HasLeaderNode(false);	//added 9-16-21
        	//nodes.push_back(p_node);
        	//node_index++;
		kk = kk + 0.5;
	}
	else	//added 9-16-21  cell mesh node accounting
	{
		Node<2>* p_node = new Node<2>(node_index, true, kk - (numelementsacross/4), ycoordinate);
		p_node->HasFocal(false);
		p_node->HasFocal1(false);
		if (i == 0 || i == numelementsacross)
		{
			p_node->HasLeaderNode(true);
			p_node->SetProtrusiveForce(true);		//added 12-7-21
			p_node->SetUpdateAdhesion1(false);		//added 12-13-21
		}
		else 
		{
			p_node->HasLeaderNode(false);
		}
        	nodes.push_back(p_node);
        	node_index++;
		kk = kk + 0.5;
		//std::cout<<"isleadernode"<<p_node->IsLeaderNode()<<std::endl;
		
	}

    }
//exit(0);

//********************************************************************************************************************************************
//************************************************CREATING ELEMENTS*****************************************************************
    
	if (isCollagen == false && isSubstrate == false)      //cell mesh condition
	{
		for (unsigned i=0; i<numElementsAcross; i++)
        	{
			node_indices[0] = i;
			node_indices[1] = i + 1;
			for (unsigned k=0; k<2; k++)
			{
				element_nodes1.push_back(nodes[node_indices[k]]);
			}
			element_index = i;
			elements.push_back(new VertexElement<2,2>(element_index, element_nodes1));
			element_nodes1.clear();						
		}
		std::cout<<"cell_elements.size()\t"<<elements.size()<<std::endl;
	}		

	else if (isCollagen == true)	 //collagen mesh, long and short fiber condition
	{
		std::cout<<"inside_collagen_condition"<<std::endl;
		//unsigned numElementsAcross_1 = nodes.size() - 1;		
		if (isLongCollagen == true)	//long fiber condition
		{
			unsigned elem_iter = 0;
			
			
			//std::cout<<"nodes.size()\t"<<nodes.size()<<std::endl;
			//exit(0);		
			
			//for (unsigned i=0; i<numElementsAcross_1-numNodesRemoved; i++)	//commented 12-6-21
			for (unsigned i=0; i<nodes.size()-1; i++)	//added 12-6-21
        		{
				//std::cout<<"i not included"<<((numElementsAcross_1 + 1 - numNodesRemoved)/2 - 1)<<std::endl;
				//if (i != ((numElementsAcross_1 + 1 - numNodesRemoved)/2 - 1))		//not registering middle element in the list of elements		//commented 12-6-21
				if (i != ((nodes.size())/2 - 1))		//added 12-6-21
				{					
					node_indices[0] = i;
					node_indices[1] = i + 1;
					//std::cout<<"node_indices[0] is\t"<<node_indices[0]<<"\t"<<"node_indices[1]\t"<<node_indices[1]<<std::endl;
					for (unsigned k=0; k<2; k++)
					{
						element_nodes1.push_back(nodes[node_indices[k]]);
					}
					//std::cout<<"element index is\t"<<i<<"\tnode_indices[0]\t"<<node_indices[0]<<"\t"<<"node_indices[1]\t"<<node_indices[1]<<std::endl;
				
					element_index = elem_iter;
					elements.push_back(new VertexElement<2,2>(element_index, element_nodes1));
					element_nodes1.clear();
					elem_iter++;
				}						
			}
		}

		//
		
		/*//commented 12-25-21
		else   //short fiber condition
		{
			//unsigned counter = 1;		//added 12-25-21
			unsigned elem_iter = 0;
			//for (unsigned i=0; i<numElementsAcross_1-numNodesRemoved; i=i+2)	//commented 12-25-21
			for (unsigned i=0; i<nodes.size()-1; i=i+2)	//added 12-25-21
        		{
				if(i!= ((nodes.size())/2 - 1))
				{
				node_indices[0] = i;
				node_indices[1] = i + 1;
				//std::cout<<"element_indices is\t"<<elem_iter<<"\t"<<"node_indices[0] is\t"<<node_indices[0]<<"\t"<<"node_indices[1]\t"<<node_indices[1]<<std::endl;
				for (unsigned k=0; k<2; k++)
				{
					element_nodes1.push_back(nodes[node_indices[k]]);
				}
				element_index = elem_iter;
				elements.push_back(new VertexElement<2,2>(element_index, element_nodes1));
				element_nodes1.clear();	
				elem_iter++;
				}					
			}
			//std::cout<<"check\t"<<std::endl;	
		}*/
	
		//added 12-25-21
		else   //short fiber condition
		{
			//unsigned counter = 1;		//added 12-25-21
			unsigned rem = (nodes.size()/2)%shortFiberConnectivity;
			unsigned elem_iter = 0;
			//std::cout<<"nodes.size())/2 - rem\t"<<nodes.size()/2<<std::endl;
			//exit(0);
			//for (unsigned i=0; i<numElementsAcross_1-numNodesRemoved; i=i+2)	//commented 12-25-21
			std::vector<unsigned> vect;
			for (unsigned i=0; i<nodes.size()/2-1; i=i+1)	//added 12-25-21
        		{
				//std::cout<<"i\t"<<i<<"\tnodes.size()\t"<<nodes.size()<<"\tnumElementsAcross_1-numNodesRemoved\t"<<numElementsAcross_1-numNodesRemoved<<std::endl;
				//std::cout<<"i\t"<<i<<std::endl;
				if(i < ((nodes.size())/2 - rem))
				{
					if((i+1)%(shortFiberConnectivity) != 0 && i != ((nodes.size())/2 - 1))
					{ 					
						node_indices[0] = i;
						node_indices[1] = i + 1;
						//std::cout<<"element_indices is\t"<<elem_iter<<"\t"<<"node_indices[0] is\t"<<node_indices[0]<<"\t"<<"node_indices[1]\t"<<node_indices[1]<<std::endl;
						for (unsigned k=0; k<2; k++)
						{
							element_nodes1.push_back(nodes[node_indices[k]]);
						}
						vect.push_back(node_indices[1]);
						element_index = elem_iter;
						elements.push_back(new VertexElement<2,2>(element_index, element_nodes1));
						element_nodes1.clear();	
						elem_iter++;
					}
				}
				
			}
			std::reverse(vect.begin(), vect.end());
			for (typename std::vector<unsigned>::iterator iter = vect.begin(); iter != vect.end(); ++iter)
			{
				node_indices[0] = nodes.size()/2 + (nodes.size()/2 - 1) - *iter;
				node_indices[1] = node_indices[0] + 1;
				for (unsigned k=0; k<2; k++)
				{
					element_nodes1.push_back(nodes[node_indices[k]]);
				}
				//std::cout<<"*iter\t"<<*iter<<"\tnodes.size()\t"<<nodes.size()<<std::endl;
				//std::cout<<"element_indices is\t"<<elem_iter<<"\t"<<"node_indices[0] is\t"<<node_indices[0]<<"\t"<<"node_indices[1]\t"<<node_indices[1]<<std::endl;
				element_index = elem_iter;
				elements.push_back(new VertexElement<2,2>(element_index, element_nodes1));
				element_nodes1.clear();	
				elem_iter++;
			}
		}
		/*
				else
				{	
					if((i+2)%(shortFiberConnectivity) != 0)
					{ 					
						node_indices[0] = i;
						node_indices[1] = i + 1;
						//std::cout<<"element_indices is\t"<<elem_iter<<"\t"<<"node_indices[0] is\t"<<node_indices[0]<<"\t"<<"node_indices[1]\t"<<node_indices[1]<<std::endl;
						for (unsigned k=0; k<2; k++)
						{
							element_nodes1.push_back(nodes[node_indices[k]]);
						}
						element_index = elem_iter;
						elements.push_back(new VertexElement<2,2>(element_index, element_nodes1));
						element_nodes1.clear();	
						elem_iter++;
					}
				}				
			}
			//std::cout<<"check\t"<<std::endl;	
		}*/
		//std::cout<<"collagen_elements.size()\t"<<elements.size()<<std::endl;
	}
	else		//substrate mesh condition
	{
		//unsigned elem_iter = 0;
		unsigned numNodes = nodes.size();
		unsigned numElementsAcross_1 = numNodes - 1;	
		for (unsigned i=0; i<numElementsAcross_1; i++)
        	{
			node_indices[0] = i;
			node_indices[1] = i + 1;
			for (unsigned k=0; k<2; k++)
			{
				element_nodes1.push_back(nodes[node_indices[k]]);
			}
			element_index = i;
			elements.push_back(new VertexElement<2,2>(element_index, element_nodes1));
			element_nodes1.clear();						
		}
		std::cout<<"substrate_elements.size()\t"<<elements.size()<<std::endl;
		
	}		
//exit(0);

    mpMesh = new MutableVertexMesh<2,2>(nodes, elements, cellRearrangementThreshold, t2Threshold);

    // Scale the mesh so that each element's area takes the value elementArea
    mpMesh->Scale(sqrt(elementArea*2.0/sqrt(3.0)), sqrt(elementArea*2.0/sqrt(3.0)));
}

StraightVertexMeshGenerator_ab2::~StraightVertexMeshGenerator_ab2()
{
    delete mpMesh;
}

MutableVertexMesh<2,2>* StraightVertexMeshGenerator_ab2::GetMesh()
{
    return mpMesh;
}
