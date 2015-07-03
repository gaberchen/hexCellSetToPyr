/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This OpenFOAM utility is distributed in the hope that it will be useful, 
    but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    hexCellSetToPyr

Description
    

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "boundaryCutter.H"
#include "cellSplitter.H"
#include "edgeCollapser.H"
#include "meshTools.H"
#include "Pair.H"
#include "globalIndex.H"
#include "cellSet.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "addOverwriteOption.H"
    argList::addOption
    (
        "set",
        "name",
        "split cells from specified cellSet only"
    );

#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createPolyMesh.H"
    const word oldInstance = mesh.pointsInstance();

    const bool overwrite = args.optionFound("overwrite");

    Info<< "Reading modifyMeshDict\n" << endl;

    bool validInputs = true;

    Info<< nl << "Looking up cells to convert to pyramids around"
        << " cell centre ..." << nl << endl;
    // Read cells in cellSet
    const bool readSet   = args.optionFound("set");
   
    // Read cells to cut from cellSet
    cellSet set(mesh, args["set"]);
    List<label> cellLabels(set.toc());
    Info<< "    Found " << cellLabels.size() << " cells in cellset " << args["set"] << endl;

    // Loop over cells in cellSet and add cell centers to topology change
    Map<point> cellToPyrCentre(cellLabels.size());


    if (readSet)
    {
        forAll(cellLabels, i)
        {                        
            // Load cell in cell set
            label cellI = cellLabels[i];

            // Insert cell center    
            if (
                    cellI == -1 ||
                    !cellToPyrCentre.insert(cellI,mesh.cellCentres()[cellI])
               )
            {
                Info<< "Could not insert mesh cell " << cellI
                    << " for input point " << mesh.cellCentres()[cellI] << nl
                    << "Perhaps the cell is already marked for splitting?" << endl;

                validInputs = false;
            }
        }
    }
    else
    {
        Info<< "cellSet must be provided using the -set option. Exiting."
            << endl;
        validInputs = false;
    }


    if (!validInputs)
    {
        Info<< nl << "There was a problem in one of the inputs in the"
            << " dictionary. Not modifying mesh." << endl;
    }
    else if (cellToPyrCentre.size())
    {
        Info<< nl << "All input cells located. Modifying mesh." << endl;

        // Mesh change engine
        cellSplitter cutter(mesh);

        // Topo change container
        polyTopoChange meshMod(mesh);

        // Insert commands into meshMod
        cutter.setRefinement(cellToPyrCentre, meshMod);

        // Do changes
        autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

        if (morphMap().hasMotionPoints())
        {
            mesh.movePoints(morphMap().preMotionPoints());
        }

        cutter.updateMesh(morphMap());

        if (!overwrite)
        {
            runTime++;
        }
        else
        {
            mesh.setInstance(oldInstance);
        }

        // Write resulting mesh
        Info<< "Writing modified mesh to time " << runTime.timeName() << endl;
        mesh.write();
    }

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
