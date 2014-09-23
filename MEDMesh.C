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

\*---------------------------------------------------------------------------*/

#include "MEDMesh.H"
#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "globalMeshData.H"
#include "PstreamCombineReduceOps.H"
#include "processorPolyPatch.H"
#include "cellModeller.H"
#include "IOmanip.H"
#include "globalIndex.H"
#include "mapDistribute.H"
#include "stringListOps.H"

#include <fstream>
#include <med.h>

// * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * * //
const Foam::label Foam::MEDMesh::foamToMEDNodeAddr[4][8] =
{
    { 5, 6, 7, 4, 1, 2, 3, 0 },     // 11 = pro-STAR hex
    { 1, 0, 2, 4, 3, 5,-1,-1 },    // 12 = pro-STAR prism
    { 1, 0, 2, 3,-1,-1,-1,-1 },   // 13 = pro-STAR tetra
    { 3, 2, 1, 0, 4,-1,-1,-1 }     // 14 = pro-STAR pyramid
};


void Foam::MEDMesh::correct()
{
    patchPartOffset_ = 2;
    meshCellSets_.setSize(mesh_.nCells());

    boundaryFaceSets_.setSize(mesh_.boundary().size());
    allPatchNames_.clear();
    patchNames_.clear();
    nPatchPrims_ = 0;
    faceZoneFaceSets_.setSize(mesh_.faceZones().size());
    faceZoneNames_.clear();
    nFaceZonePrims_ = 0;
    boundaryFaceToBeIncluded_.clear();

    if (!noPatches_)
    {
        // Patches are output. Check that they're synced.
        mesh_.boundaryMesh().checkParallelSync(true);

        allPatchNames_ = mesh_.boundaryMesh().names();
        if (Pstream::parRun())
        {
            allPatchNames_.setSize
            (
                mesh_.boundary().size()
              - mesh_.globalData().processorPatches().size()
            );
        }

        if (patches_)
        {
            if (patchPatterns_.empty())
            {
                forAll(allPatchNames_, nameI)
                {
                    patchNames_.insert(allPatchNames_[nameI]);
                }
            }
            else
            {
                // Find patch names which match that requested at command-line
                forAll(allPatchNames_, nameI)
                {
                    const word& patchName = allPatchNames_[nameI];
                    if (findStrings(patchPatterns_, patchName))
                    {
                        patchNames_.insert(patchName);
                    }
                }
            }
        }
    }

    if (patchNames_.size())
    {
        // no internalMesh
        patchPartOffset_ = 1;
    }
    else
    {
        const cellShapeList& cellShapes = mesh_.cellShapes();

        const cellModel& tet = *(cellModeller::lookup("tet"));
        const cellModel& pyr = *(cellModeller::lookup("pyr"));
        const cellModel& prism = *(cellModeller::lookup("prism"));
        const cellModel& wedge = *(cellModeller::lookup("wedge"));
        const cellModel& hex = *(cellModeller::lookup("hex"));



        // Count the shapes
        labelList& tets = meshCellSets_.tets;
        labelList& pyrs = meshCellSets_.pyrs;
        labelList& prisms = meshCellSets_.prisms;
        labelList& wedges = meshCellSets_.wedges;
        labelList& hexes = meshCellSets_.hexes;
        labelList& polys = meshCellSets_.polys;

        label nTets = 0;
        label nPyrs = 0;
        label nPrisms = 0;
        label nWedges = 0;
        label nHexes = 0;
        label nPolys = 0;

        forAll(cellShapes, cellI)
        {
            const cellShape& cellShape = cellShapes[cellI];
            const cellModel& cellModel = cellShape.model();

            if (cellModel == tet)
            {
                tets[nTets++] = cellI;
            }
            else if (cellModel == pyr)
            {
                pyrs[nPyrs++] = cellI;
            }
            else if (cellModel == prism)
            {
                prisms[nPrisms++] = cellI;
            }
            else if (cellModel == wedge)
            {
                wedges[nWedges++] = cellI;
            }
            else if (cellModel == hex)
            {
                hexes[nHexes++] = cellI;
            }
            else
            {
                polys[nPolys++] = cellI;
            }
        }

        tets.setSize(nTets);
        pyrs.setSize(nPyrs);
        prisms.setSize(nPrisms);
        wedges.setSize(nWedges);
        hexes.setSize(nHexes);
        polys.setSize(nPolys);

        meshCellSets_.nTets = nTets;
        reduce(meshCellSets_.nTets, sumOp<label>());

        meshCellSets_.nPyrs = nPyrs;
        reduce(meshCellSets_.nPyrs, sumOp<label>());

        meshCellSets_.nPrisms = nPrisms;
        reduce(meshCellSets_.nPrisms, sumOp<label>());

        meshCellSets_.nHexesWedges = nWedges+nHexes;
        reduce(meshCellSets_.nHexesWedges, sumOp<label>());

        meshCellSets_.nPolys = nPolys;
        reduce(meshCellSets_.nPolys, sumOp<label>());


        // Determine parallel shared points
        globalPointsPtr_ = mesh_.globalData().mergePoints
        (
            pointToGlobal_,
            uniquePointMap_
        );
    }

    if (!noPatches_)
    {
        forAll(mesh_.boundary(), patchi)
        {
            if (mesh_.boundary()[patchi].size())
            {
                const polyPatch& p = mesh_.boundaryMesh()[patchi];

                labelList& tris = boundaryFaceSets_[patchi].tris;
                labelList& quads = boundaryFaceSets_[patchi].quads;
                labelList& polys = boundaryFaceSets_[patchi].polys;

                tris.setSize(p.size());
                quads.setSize(p.size());
                polys.setSize(p.size());

                label nTris = 0;
                label nQuads = 0;
                label nPolys = 0;

                forAll(p, faceI)
                {
                    const face& f = p[faceI];

                    if (f.size() == 3)
                    {
                        tris[nTris++] = faceI;
                    }
                    else if (f.size() == 4)
                    {
                        quads[nQuads++] = faceI;
                    }
                    else
                    {
                        polys[nPolys++] = faceI;
                    }
                }

                tris.setSize(nTris);
                quads.setSize(nQuads);
                polys.setSize(nPolys);
            }
        }
    }

    forAll(allPatchNames_, patchi)
    {
        const word& patchName = allPatchNames_[patchi];
        nFacePrimitives nfp;

        if (patchNames_.empty() || patchNames_.found(patchName))
        {
            if (mesh_.boundary()[patchi].size())
            {
                nfp.nTris   = boundaryFaceSets_[patchi].tris.size();
                nfp.nQuads  = boundaryFaceSets_[patchi].quads.size();
                nfp.nPolys  = boundaryFaceSets_[patchi].polys.size();
            }
        }

        reduce(nfp.nTris, sumOp<label>());
        reduce(nfp.nQuads, sumOp<label>());
        reduce(nfp.nPolys, sumOp<label>());

        nPatchPrims_.insert(patchName, nfp);
    }

    nTris_ = 0;
    nQuads_ = 0;
    nPolys_ = 0;
    forAll(allPatchNames_, patchi)
    {
        const word& patchName = allPatchNames_[patchi];
        nTris_ += nPatchPrims_[patchName].nTris;
        nQuads_ += nPatchPrims_[patchName].nQuads;
        nPolys_ += nPatchPrims_[patchName].nPolys;
    }

    // faceZones
    if (faceZones_)
    {
        wordList faceZoneNamesAll = mesh_.faceZones().names();
        // Need to sort the list of all face zones since the index may vary
        // from processor to processor...
        sort(faceZoneNamesAll);

        // Find faceZone names which match that requested at command-line
        forAll(faceZoneNamesAll, nameI)
        {
            const word& zoneName = faceZoneNamesAll[nameI];
            if (findStrings(faceZonePatterns_, zoneName))
            {
                faceZoneNames_.insert(zoneName);
            }
        }

        // Build list of boundary faces to be exported
        boundaryFaceToBeIncluded_.setSize
        (
            mesh_.nFaces()
          - mesh_.nInternalFaces(),
            1
        );

        forAll(mesh_.boundaryMesh(), patchI)
        {
            const polyPatch& pp = mesh_.boundaryMesh()[patchI];
            if
            (
                isA<processorPolyPatch>(pp)
             && !refCast<const processorPolyPatch>(pp).owner()
            )
            {
                label bFaceI = pp.start()-mesh_.nInternalFaces();
                forAll(pp, i)
                {
                    boundaryFaceToBeIncluded_[bFaceI++] = 0;
                }
            }
        }

        // Count face types in each faceZone
        forAll(faceZoneNamesAll, zoneI)
        {
            const word& zoneName = faceZoneNamesAll[zoneI];
            const label faceZoneId = mesh_.faceZones().findZoneID(zoneName);

            const faceZone& fz = mesh_.faceZones()[faceZoneId];

            if (fz.size())
            {
                labelList& tris = faceZoneFaceSets_[faceZoneId].tris;
                labelList& quads = faceZoneFaceSets_[faceZoneId].quads;
                labelList& polys = faceZoneFaceSets_[faceZoneId].polys;

                tris.setSize(fz.size());
                quads.setSize(fz.size());
                polys.setSize(fz.size());

                label nTris = 0;
                label nQuads = 0;
                label nPolys = 0;

                label faceCounter = 0;

                forAll(fz, i)
                {
                    label faceI = fz[i];

                    // Avoid counting faces on processor boundaries twice
                    if (faceToBeIncluded(faceI))
                    {
                        const face& f = mesh_.faces()[faceI];

                        if (f.size() == 3)
                        {
                            tris[nTris++] = faceCounter;
                        }
                        else if (f.size() == 4)
                        {
                            quads[nQuads++] = faceCounter;
                        }
                        else
                        {
                            polys[nPolys++] = faceCounter;
                        }

                        ++faceCounter;
                    }
                }

                tris.setSize(nTris);
                quads.setSize(nQuads);
                polys.setSize(nPolys);
            }
        }

        forAll(faceZoneNamesAll, zoneI)
        {
            const word& zoneName = faceZoneNamesAll[zoneI];
            nFacePrimitives nfp;
            const label faceZoneId = mesh_.faceZones().findZoneID(zoneName);

            if (faceZoneNames_.found(zoneName))
            {
                if
                (
                    faceZoneFaceSets_[faceZoneId].tris.size()
                 || faceZoneFaceSets_[faceZoneId].quads.size()
                 || faceZoneFaceSets_[faceZoneId].polys.size()
                )
                {
                    nfp.nTris   = faceZoneFaceSets_[faceZoneId].tris.size();
                    nfp.nQuads  = faceZoneFaceSets_[faceZoneId].quads.size();
                    nfp.nPolys  = faceZoneFaceSets_[faceZoneId].polys.size();
                }
            }

            reduce(nfp.nTris, sumOp<label>());
            reduce(nfp.nQuads, sumOp<label>());
            reduce(nfp.nPolys, sumOp<label>());

            nFaceZonePrims_.insert(zoneName, nfp);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MEDMesh::MEDMesh
(
    const fvMesh& mesh,
    const bool noPatches,

    const bool patches,
    const wordReList& patchPatterns,

    const bool faceZones,
    const wordReList& faceZonePatterns,

    const bool binary
)
:
    mesh_(mesh),
    noPatches_(noPatches),
    patches_(patches),
    patchPatterns_(patchPatterns),
    faceZones_(faceZones),
    faceZonePatterns_(faceZonePatterns),
    binary_(binary),
    meshCellSets_(mesh.nCells())
{
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::MEDMesh::~MEDMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::MEDMesh::faceToBeIncluded(const label faceI) const
{
    bool res = false;

    if (mesh_.isInternalFace(faceI))
    {
        res = true;
    }
    else
    {
        res = boundaryFaceToBeIncluded_[faceI-mesh_.nInternalFaces()];
    }

    return res;
}


void Foam::MEDMesh::barrier()
{
    label appI = 0;
    reduce(appI,maxOp<label>());
}


Foam::cellShapeList Foam::MEDMesh::map
(
    const cellShapeList& cellShapes,
    const labelList& prims,
    const labelList& pointToGlobal
) const
{
    cellShapeList mcsl(prims.size());

    forAll(prims, i)
    {
        mcsl[i] = cellShapes[prims[i]];
        inplaceRenumber(pointToGlobal, mcsl[i]);
    }

    return mcsl;
}


Foam::cellShapeList Foam::MEDMesh::map
(
    const cellShapeList& cellShapes,
    const labelList& hexes,
    const labelList& wedges,
    const labelList& pointToGlobal
) const
{
    cellShapeList mcsl(hexes.size() + wedges.size());

    forAll(hexes, i)
    {
        mcsl[i] = cellShapes[hexes[i]];
        inplaceRenumber(pointToGlobal, mcsl[i]);
    }

    label offset = hexes.size();

    const cellModel& hex = *(cellModeller::lookup("hex"));
    labelList hexLabels(8);

    forAll(wedges, i)
    {
        const cellShape& cellPoints = cellShapes[wedges[i]];

        hexLabels[0] = cellPoints[0];
        hexLabels[1] = cellPoints[1];
        hexLabels[2] = cellPoints[0];
        hexLabels[3] = cellPoints[2];
        hexLabels[4] = cellPoints[3];
        hexLabels[5] = cellPoints[4];
        hexLabels[6] = cellPoints[6];
        hexLabels[7] = cellPoints[5];

        mcsl[i + offset] = cellShape(hex, hexLabels);
        inplaceRenumber(pointToGlobal, mcsl[i + offset]);
    }

    return mcsl;
}


void Foam::MEDMesh::writePrims
(
    const med_geometry_type key,
    const cellShapeList& cellShapes,
    const label *mapping,
    const char *meshname,
    const int medfileid
) const
{
    // Create a temp int array
    if (cellShapes.size())
    {
        // All the cellShapes have the same number of elements!
        int nelems = cellShapes.size();
        int nnodes = cellShapes[0].size();
        int *temp = new med_int[nelems*nnodes];
        
        forAll(cellShapes, i)
        {
            const cellShape& cellPoints = cellShapes[i];
            forAll(cellPoints, pointI)
            {
                temp[i*nnodes+mapping[pointI]] = cellPoints[pointI] + 1;
            }
        }
        med_err ret = MEDmeshElementConnectivityWr(
            medfileid, meshname,MED_NO_DT,MED_NO_IT,MED_NO_DT, MED_CELL,
            key, MED_NODAL,MED_FULL_INTERLACE, nelems,temp);

        delete[] temp;
        if (ret == 0)
            Info << "wrote " << nelems << " prisms with each "<< nnodes <<" nodes" << endl;
        else
            Info << "failed to write prisms with med geotype "<<key<<endl;

    }
}

void Foam::MEDMesh::writeAllPolyhedrons
(
    const labelList& pointToGlobal,
    const char* meshname,
    const int medfileid
) const
{
    if (meshCellSets_.nPolys)
    {
        const cellList& cellFaces = mesh_.cells();
        const labelList& faceOwner = mesh_.faceOwner();

        // Renumber faces to use global point numbers
        faceList faces(mesh_.faces());
        forAll(faces, i)
        {
            inplaceRenumber(pointToGlobal, faces[i]);
        }
        // Number of faces for each poly cell
        const labelList& polys = meshCellSets_.polys;
        const int nelems = polys.size();
        med_int *temp_ifdx = new med_int[nelems+1];
        med_int c_ifdx = 1;
        temp_ifdx[0] = c_ifdx;
        forAll(polys,i)
        {
            c_ifdx += cellFaces[polys[i]].size();
            temp_ifdx[i+1] = c_ifdx;
        }
        // Number of points for each face of the above list
        med_int *temp_idx = new med_int[c_ifdx];
        med_int c_idx = 1;
        int c = 0;
        temp_idx[c++] = c_idx;
        forAll(polys,i)
        {
            const labelList& cf = cellFaces[polys[i]];
            forAll(cf, faceI)
            {
                c_idx += faces[cf[faceI]].size();
                temp_idx[c++] = c_idx;
            }
        }
        // List of points id for each face of the above list
        med_int *temp_conns = new med_int[c_idx];
        c = 0;
        forAll(polys,i)
        {
            const labelList& cf = cellFaces[polys[i]];
            forAll(cf, faceI)
            {
                const label faceId = cf[faceI];
                const face& f = faces[faceId];  // points of face (in global points)
                const label np = f.size();
                bool reverseOrder = false;
                if (faceId >= faceOwner.size())
                {
                    // Boundary face.
                    // Nothing should be done for processor boundary.
                    // The current cell always owns them. Given that we
                    // are reverting the
                    // order when the cell is the neighbour to the face,
                    // the orientation of
                    // all the boundaries, no matter if they are "real"
                    // or processorBoundaries, is consistent.
                }
                else
                {
                    if (faceOwner[faceId] != polys[i])
                        reverseOrder = true;
                }
                // If the face owner is the current cell, write the points
                // in the standard order.
                // If the face owner is not the current cell, write the points
                // in reverse order.
                // MED prefers to have all the faces of an nfaced cell
                // oriented in the same way.
                forAll(f, pointI)
                {
                    if (reverseOrder)
                        temp_conns[c++] = f[np-1-pointI] + 1;
                    else
                        temp_conns[c++] = f[pointI] + 1;
                }
            }
        }

        // now write the data to the file
        med_err ret = MEDmeshPolyhedronWr(
                medfileid, meshname,MED_NO_DT,MED_NO_IT,MED_NO_DT, MED_CELL,
                MED_NODAL,nelems+1,temp_ifdx,c_ifdx,temp_idx,temp_conns);
        if (ret!=0)
            Info << "an error occured when writing the polyhedrons" << endl;
        else
            Info << "wrote " << nelems << " polyhedrons" << endl;
        delete[] temp_ifdx;
        delete[] temp_idx;
        delete[] temp_conns;

    }
}


void Foam::MEDMesh::writeAllPrims
(
    const med_geometry_type key,
    const int nPrisms,
    const cellShapeList& cellShapes,
    const label *mapping,
    const char* meshname,
    const int medfileid
) const
{
    if (nPrisms)
    {
        writePrims(key,cellShapes,mapping, meshname,medfileid);

    }
}


void Foam::MEDMesh::writeFamilies
(
 const med_geometry_type key,
 const int nelems,
 const med_int *families,
 const char* meshname,
 const int medfileid
 )const
{
    if (nelems > 0)
    {
        med_err ret = MEDmeshEntityFamilyNumberWr(
                medfileid, meshname,MED_NO_DT,MED_NO_IT, MED_CELL,
                key,nelems,families);
        if (ret!=0)
            Info << "failed to write families for " << key << endl;
        else
            Info << "wrote families for " << key << endl;
    }
}

void Foam::MEDMesh::writeFacePrims
(
    const med_geometry_type key,
    const faceList& patchFaces,
    const char* meshname,
    const int medfileid
) const
{
    if (patchFaces.size())
    {
        int nelems = patchFaces.size();
        int nnodes = patchFaces[0].size();
        med_int *temp = new med_int[nelems*nnodes];
        int c = 0;
        forAll(patchFaces, i)
        {
            const face& patchFace = patchFaces[i];
            forAll(patchFace, pointI)
            {
                temp[c++] = patchFace[pointI] + 1;
            }
        }
        med_err ret = MEDmeshElementConnectivityWr(
                medfileid, meshname,MED_NO_DT,MED_NO_IT,MED_NO_DT, MED_CELL,
                key, MED_NODAL,MED_FULL_INTERLACE, nelems,temp);

        delete[] temp;
        if (ret == 0)
            Info << "wrote " << nelems << " face prisms with each "<< nnodes <<" nodes" << endl;
        else
            Info << "failed to write prisms with med geotype "<<key<<endl;
    }


}


void Foam::MEDMesh::writeAllFacePrims
(
    const med_geometry_type key,
    const labelList& prims,
    const label nPrims,
    const faceList& patchFaces,
    const char* meshname,
    const int medfileid
) const
{
    if (nPrims)
    {
        writeFacePrims(
                key,
                UIndirectList<face>(patchFaces, prims)(),
                meshname,
                medfileid);

    }
}

void Foam::MEDMesh::writeAllPolygons
(
    const labelList& prims,
    const label nPrims,
    const faceList& patchFaces,
    const char* meshname,
    const int medfileid
) const
{
    if (nPrims)
    {
        int nelems = patchFaces.size();
        med_int *temp_idx = new med_int[nelems+1];
        int c_idx = 1;
        temp_idx[0] = c_idx;
        forAll(patchFaces, i)
        {
            c_idx += patchFaces[i].size();
            temp_idx[i+1] = c_idx;
        }
        med_int *temp_conns = new med_int[c_idx+1];
        int c = 0;
        forAll(patchFaces, i)
        {
            const face& patchFace = patchFaces[i];
            forAll(patchFace, pointI)
                temp_conns[c++] = patchFace[pointI] + 1;
            
        }

        med_err ret = MEDmeshPolygonWr(
                medfileid,meshname,MED_NO_DT,MED_NO_IT,MED_NO_DT,MED_CELL,
                MED_NODAL,nelems+1,temp_idx,temp_conns);

        if (ret!=0)
            Info << "an error occured when writing the polygons" << endl;
        else
            Info << "wrote " << nelems << " polygons" << endl;
        delete[] temp_idx;
        delete[] temp_conns;
    }
}


void Foam::MEDMesh::writeAllPoints
(
    const pointField& uniquePoints,
    const label nPoints,
    const char *meshname,
    const int medfileid
) const
{
    barrier();

    // create a temporary array which is continuos 
    med_float * tmp = new med_float[nPoints*3*sizeof(med_float)];

    forAll(uniquePoints, pointI)
    {
        tmp[pointI*3] = uniquePoints[pointI].x();
        tmp[pointI*3+1] = uniquePoints[pointI].y();
        tmp[pointI*3+2] = uniquePoints[pointI].z();
    }

    MEDmeshNodeCoordinateWr(
            medfileid,meshname,MED_NO_DT,MED_NO_IT,MED_NO_DT,
            MED_FULL_INTERLACE,nPoints,tmp);
    // delete the temporary array
    delete[] tmp;

}


void Foam::MEDMesh::write
(
    const label timeIndex,
    const bool meshMoving,
    const char * meshname,
    const int medfileid
) const
{
    const Time& runTime = mesh_.time();
    const cellShapeList& cellShapes = mesh_.cellShapes();


   if (patchNames_.empty())
    {
        label nPoints = globalPoints().size();

        const pointField uniquePoints(mesh_.points(), uniquePointMap_);

        writeAllPoints
        (
            uniquePoints,
            nPoints,
            meshname,
            medfileid
        );

        writeAllPrims
        (
            MED_HEXA8,
            meshCellSets_.nHexesWedges,
            map         // Rewrite cellShapes to global numbering
            (
                cellShapes,
                meshCellSets_.hexes,
                meshCellSets_.wedges,
                pointToGlobal_
            ),
            foamToMEDNodeAddr[0],
            meshname,
            medfileid
        );

        writeAllPrims
        (
            MED_PENTA6,
            meshCellSets_.nPrisms,
            map(cellShapes, meshCellSets_.prisms, pointToGlobal_),
            foamToMEDNodeAddr[1],
            meshname,
            medfileid
        );

        writeAllPrims
        (
            MED_PYRA5,
            meshCellSets_.nPyrs,
            map(cellShapes, meshCellSets_.pyrs, pointToGlobal_),
            foamToMEDNodeAddr[3],
            meshname,
            medfileid
        );

        writeAllPrims
        (
            MED_TETRA4,
            meshCellSets_.nTets,
            map(cellShapes, meshCellSets_.tets, pointToGlobal_),
            foamToMEDNodeAddr[2],
            meshname,
            medfileid
        );

        writeAllPolyhedrons
        (
            pointToGlobal_,
            meshname,
            medfileid
        );
    }

    // now write the patches
    // first we collect all the different facestypes oft the patches and write
    // them to the med file
    // then we write the patchnumbers as familynumbers
    // and finally we write the familyinfo with the patchname as the groupname
    labelList *trisPtr = new labelList(nTris_);
    labelList &tris = *trisPtr;
    med_int *trifams = new med_int[nTris_];
    int tric = 0;

    labelList *quadsPtr = new labelList(nQuads_);
    labelList &quads = *quadsPtr;
    med_int *quadfams = new med_int[nQuads_];
    int quadc = 0;

    labelList *polysPtr = new labelList(nPolys_);
    labelList &polys = *polysPtr;
    med_int *polyfams = new med_int[nPolys_];
    int polyc = 0;

    const faceList& faces = mesh_.faces();

    forAll(allPatchNames_, patchi)
    {
        const word& patchName = allPatchNames_[patchi];
        const polyPatch& p = mesh_.boundaryMesh()[patchi];
        const int end = p.start()+p.size();
        for(int i=p.start(); i < end ; i++)
        {
            if (faces[i].size() == 3)
            {
                tris[tric] = i;
                trifams[tric++] = -1*(patchi+1);
            }
            else if (faces[i].size() == 4)
            {
                quads[quadc] = i;
                quadfams[quadc++] = -1*(patchi+1);
            }
            else
            {
                polys[polyc] = i;
                polyfams[polyc++] = -1*(patchi+1);
            }
        }

    }

    writeAllFacePrims(MED_TRIA3,tris,tric,faces,meshname,medfileid);
    writeAllFacePrims(MED_QUAD4,quads,quadc,faces,meshname,medfileid);
    writeAllFacePrims(MED_POLYGON,polys,polyc,faces,meshname,medfileid);
    writeFamilies(MED_TRIA3,tric,trifams,meshname,medfileid);
    writeFamilies(MED_QUAD4,quadc,quadfams,meshname,medfileid);
    writeFamilies(MED_POLYGON,polyc,polyfams,meshname,medfileid);

    med_err ret;
    int famnum;
    char familyname[MED_NAME_SIZE];
    char patchname[MED_LNAME_SIZE];
    forAll(allPatchNames_, patchi)
    {
        const word& patchName = allPatchNames_[patchi];
        famnum = -1*(patchi+1);
        std::string fname = "FAM_"+std::to_string(famnum);
        for(unsigned int i=0;i<MED_NAME_SIZE;i++)
        {
            if (i<fname.size())
                familyname[i] = fname[i];
            else
                familyname[i] = '\0';
        }
        for(unsigned int i=0;i<MED_LNAME_SIZE;i++)
        {
            if (i<patchName.size())
                patchname[i] = patchName[i];
            else
                patchname[i] = '\0';
        }
        ret = MEDfamilyCr
            (
            medfileid,
            meshname,
            familyname,
            famnum,
            1,
            patchname
            );
        if (ret != 0)
            Info << "failed to write familyinfo for patch " << patchName << endl;
        else
            Info << "wrote familyinfo for patch " << patchName << endl;
   }


  /*  // write faceZones, if requested*/
    //forAllConstIter(wordHashSet, faceZoneNames_, iter)
    //{
        //const word& faceZoneName = iter.key();

        //label faceID = mesh_.faceZones().findZoneID(faceZoneName);

        //const faceZone& fz = mesh_.faceZones()[faceID];

        //const nFacePrimitives& nfp = nFaceZonePrims_[faceZoneName];

        //if (nfp.nTris || nfp.nQuads || nfp.nPolys)
        //{
            //const labelList& tris = faceZoneFaceSets_[faceID].tris;
            //const labelList& quads = faceZoneFaceSets_[faceID].quads;
            //const labelList& polys = faceZoneFaceSets_[faceID].polys;

            //// Renumber the faceZone points/faces into unique points
            //labelList pointToGlobal;
            //labelList uniqueMeshPointLabels;
            //autoPtr<globalIndex> globalPointsPtr =
                //mesh_.globalData().mergePoints
                //(
                    //fz().meshPoints(),
                    //fz().meshPointMap(),
                    //pointToGlobal,
                    //uniqueMeshPointLabels
                //);

            //pointField uniquePoints(mesh_.points(), uniqueMeshPointLabels);

            //// Find the list of master faces belonging to the faceZone,
            //// in local numbering
            //faceList faceZoneFaces(fz().localFaces());

            //// Count how many master faces belong to the faceZone. Is there
            //// a better way of doing this?
            //label nMasterFaces = 0;

            //forAll(fz, faceI)
            //{
                //if (faceToBeIncluded(fz[faceI]))
                //{
                    //++nMasterFaces;
                //}
            //}

            //// Create the faceList for the master faces only and fill it.
            //faceList faceZoneMasterFaces(nMasterFaces);

            //label currentFace = 0;

            //forAll(fz, faceI)
            //{
                //if (faceToBeIncluded(fz[faceI]))
                //{
                    //faceZoneMasterFaces[currentFace] = faceZoneFaces[faceI];
                    //++currentFace;
                //}
            //}

            //// Renumber the faceZone master faces
            //forAll(faceZoneMasterFaces, i)
            //{
                //inplaceRenumber(pointToGlobal, faceZoneMasterFaces[i]);
            //}

            //writeAllFacePrims
            //(
                //MED_TRIA3,
                //tris,
                //nfp.nTris,
                //faceZoneMasterFaces,
                //meshname,
                //medfileid
            //);

            //writeAllFacePrims
            //(
                //MED_QUAD4,
                //quads,
                //nfp.nQuads,
                //faceZoneMasterFaces,
                //meshname,
                //medfileid
            //);

            //writeAllPolygons
            //(
                //polys,
                //nfp.nPolys,
                //faceZoneMasterFaces,
                //meshname,
                //medfileid
            //);
        //}
    //}

}


// ************************************************************************* //
