/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    foamToMED

Description
    Translates OpenFOAM data to MED format.

    An MED part is created for the internalMesh and for each patch.

Usage
    - foamToMED [OPTION] \n
    Translates OpenFOAM data to MED format

    \param -ascii \n
    Write MED data in ASCII format instead of "C Binary"

    \param -patches patchList \n
    Specify particular patches to write.
    Specifying an empty list suppresses writing the internalMesh.

    \param -noPatches \n
    Suppress writing any patches.

    \param -faceZones zoneList \n
    Specify faceZones to write, with wildcards

    \param -cellZone zoneName \n
    Specify single cellZone to write (not lagrangian)

Note
    Parallel support for cloud data is not supported
    - writes to \a MED directory to avoid collisions with foamToMEDParts

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "IOobjectList.H"
#include "IOmanip.H"
#include "OFstream.H"

#include "MEDMesh.H"
#include "globalIndex.H"

#include "fvc.H"

#include "cellSet.H"
#include "fvMeshSubset.H"
#include "itoa.H"
#include <med.h>
#include <stdlib.h>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool inFileNameList
(
    const fileNameList& nameList,
    const word& name
)
{
    forAll(nameList, i)
    {
        if (nameList[i] == name)
        {
            return true;
        }
    }

    return false;
}



int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::addOption
    (
        "scale",
        "double",
        "specify scaling factor"
    );

    #include "setRootCase.H"

    // Check options
    const bool is_scale = args.optionFound("scale");
    double scale;
    if (is_scale)
        scale = atof(args.option("scale").c_str());
    else
       scale = 1.0; 

    #include "createTime.H"

    instantList Times = timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"

    // Start of case file header output
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const word prepend = args.globalCaseName() + '.';

    med_err ret;

    Info<< "Createing MEDfile" << endl;
    const med_idt medfile = MEDfileOpen((prepend+"med").c_str(),MED_ACC_CREAT);

    if (medfile < 0)
        Info << "failed to open the med file" << endl;

    Info<< "Writing Comment" << endl;
    std::string comment ("FOAM to MED Mesh File",MED_COMMENT_SIZE);
    ret = MEDfileCommentWr(medfile,comment.c_str());

    if (ret != 0)
        Info << "failed to write comments" << endl;

    Info<< "Creating Mesh" << endl;
    const char *meshname = args.globalCaseName().c_str();
    std::string description ("created by OpenFOAM",MED_COMMENT_SIZE);
    std::string dtunit ("s",MED_SNAME_SIZE);
    std::string axis_x ("x",MED_SNAME_SIZE);
    std::string axis_y ("y",MED_SNAME_SIZE);
    std::string axis_z ("z",MED_SNAME_SIZE);
    std::string unit_x ("m",MED_SNAME_SIZE);
    std::string unit_y ("m",MED_SNAME_SIZE);
    std::string unit_z ("m",MED_SNAME_SIZE);

    // create a mesh
    ret = MEDmeshCr
        (
            medfile,
            meshname,
            3,
            3,
            MED_UNSTRUCTURED_MESH,
            description.c_str(),
            dtunit.c_str(),
            MED_SORT_DTIT,
            MED_CARTESIAN,
            (axis_x+axis_y+axis_z).c_str(),
            (unit_x+unit_y+unit_z).c_str()
        );
    
    if (ret != 0)
        Info << "failed to create the mesh" << endl;
    //
    // Construct the MED mesh
    const bool selectedPatches = args.optionFound("patches");
    wordReList patchPatterns;
    if (selectedPatches)
    {
        patchPatterns = wordReList(args.optionLookup("patches")());
    }
    const bool selectedZones = args.optionFound("faceZones");
    wordReList zonePatterns;
    if (selectedZones)
    {
        zonePatterns = wordReList(args.optionLookup("faceZones")());
    }

    word cellZoneName;
    const bool doCellZone = args.optionReadIfPresent("cellZone", cellZoneName);

    fvMeshSubset meshSubsetter(mesh);
    if (doCellZone)
    {
        Info<< "Converting cellZone " << cellZoneName
            << " only (puts outside faces into patch "
            << mesh.boundaryMesh()[0].name()
            << ")" << endl;
        const cellZone& cz = mesh.cellZones()[cellZoneName];
        cellSet c0(mesh, "c0", labelHashSet(cz));
        meshSubsetter.setLargeCellSubset(c0, 0);
    }

    MEDMesh eMesh
    (
        (
            meshSubsetter.hasSubMesh()
          ? meshSubsetter.subMesh()
          : meshSubsetter.baseMesh()
        ),
        args.optionFound("noPatches"),
        selectedPatches,
        patchPatterns,
        selectedZones,
        zonePatterns,
        false
    );

    // Set Time to the last time before looking for the lagrangian objects
    runTime.setTime(Times.last(), Times.size()-1);

    IOobjectList objects(mesh, runTime.timeName());

    #include "checkMeshMoving.H"

    if (meshMoving)
    {
        Info<< "Detected a moving mesh (multiple polyMesh/points files)."
            << " Writing meshes for every timestep." << endl;
    }



    label nTimeSteps = 0;
    forAll(Times, timeIndex)
    {
        nTimeSteps++;
        runTime.setTime(Times[timeIndex], timeIndex);

        word timeName = itoa(timeIndex);
        word timeFile = prepend + timeName;

        Info<< "Translating time = " << runTime.timeName() << nl;

        polyMesh::readUpdateState meshState = mesh.readUpdate();
        if (timeIndex != 0 && meshSubsetter.hasSubMesh())
        {
            Info<< "Converting cellZone " << cellZoneName
                << " only (puts outside faces into patch "
                << mesh.boundaryMesh()[0].name()
                << ")" << endl;
            const cellZone& cz = mesh.cellZones()[cellZoneName];
            cellSet c0(mesh, "c0", labelHashSet(cz));
            meshSubsetter.setLargeCellSubset(c0, 0);
        }


        if (meshState != polyMesh::UNCHANGED)
        {
            eMesh.correct();
        }

        if (timeIndex == 0 || meshMoving)
        {
            eMesh.write
            (
                timeIndex,
                meshMoving,
                scale,
                meshname,
                medfile
            );
        }


        
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
