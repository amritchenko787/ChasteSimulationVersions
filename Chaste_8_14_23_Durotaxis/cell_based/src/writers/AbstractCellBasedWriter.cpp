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
#include "AbstractCellBasedWriter.hpp"
#include "SimulationTime.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCellBasedWriter<ELEMENT_DIM, SPACE_DIM>::AbstractCellBasedWriter(const std::string& rFileName)
    : mFileName(rFileName)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCellBasedWriter<ELEMENT_DIM, SPACE_DIM>::~AbstractCellBasedWriter()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedWriter<ELEMENT_DIM, SPACE_DIM>::CloseFile()
{
    mpOutStream->close();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedWriter<ELEMENT_DIM, SPACE_DIM>::OpenOutputFile(OutputFileHandler& rOutputFileHandler)
{
    mpOutStream = rOutputFileHandler.OpenOutputFile(mFileName);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedWriter<ELEMENT_DIM, SPACE_DIM>::OpenOutputFileForAppend(OutputFileHandler& rOutputFileHandler)
{
    mpOutStream = rOutputFileHandler.OpenOutputFile(mFileName, std::ios::app);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedWriter<ELEMENT_DIM, SPACE_DIM>::WriteTimeStamp()
{
    *mpOutStream << SimulationTime::Instance()->GetTime() << "\t";
	//std::cout<<"SimulationTime::Instance()->GetTime()"<<SimulationTime::Instance()->GetTime()<<std::endl;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedWriter<ELEMENT_DIM, SPACE_DIM>::WriteNewline()
{
    *mpOutStream << "\n";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedWriter<ELEMENT_DIM, SPACE_DIM>::SetFileName(std::string fileName)
{
    mFileName = fileName;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::string AbstractCellBasedWriter<ELEMENT_DIM, SPACE_DIM>::GetFileName()
{
    return mFileName;
}

// Explicit instantiation
template class AbstractCellBasedWriter<1,1>;
template class AbstractCellBasedWriter<1,2>;
template class AbstractCellBasedWriter<2,2>;
template class AbstractCellBasedWriter<1,3>;
template class AbstractCellBasedWriter<2,3>;
template class AbstractCellBasedWriter<3,3>;