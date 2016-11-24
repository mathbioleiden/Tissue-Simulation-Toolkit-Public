/* 

Copyright 1996-2006 Roeland Merks, Paulien Hogeweg

This file is part of Tissue Simulation Toolkit.

Tissue Simulation Toolkit is free software; you can redistribute
it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

Tissue Simulation Toolkit is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Tissue Simulation Toolkit; if not, write to the Free
Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
02110-1301 USA

*/
#include <vector>
#include <list>
#include <algorithm>
#include <fstream>
#include <string.h>
#include <errno.h>
#include <math.h>
#include "dish.h"
#include "sticky.h"
#include "parameter.h"
#include "info.h"
#include "crash.h"
#include "pde.h"

#include <MultiCellDS.hpp>
#include <MultiCellDS-pimpl.hpp>
#include <MultiCellDS-simpl.hpp>
//#include "pugixml.hpp"

#define EXTERNAL_OFF

extern Parameter par;

using namespace std;

Dish::Dish(void) {

  ConstructorBody();
  
  CPM=new CellularPotts(&cell, par.sizex, par.sizey);
  if (par.n_chem)
    PDEfield=new PDE(par.n_chem,par.sizex, par.sizey);
  
  // Initial cell distribution is defined by user in INIT {} block
  Init();
    
  if (par.target_area>0)
    for (std::vector<Cell>::iterator c=cell.begin();c!=cell.end();c++) {
      c->SetTargetArea(par.target_area);
    } 
  
  
}



Dish::~Dish() {
    cell.clear();
    
    delete CPM;
	
 }

void Dish::Plot(Graphics *g) {
    if (CPM)
      CPM->Plot(g);
 }


void Dish::ConstructorBody() {
  
  Cell::maxsigma=0;
  
  // Allocate the first "cell": this is the medium (tau=0)
  cell.push_back(*(new Cell(*this,0)));
  
  // indicate that the first cell is the medium
  cell.front().sigma=0; 
  cell.front().tau=0;
  
  CPM=0;
  PDEfield=0;

}


bool Dish::CellLonelyP(const Cell &c, int **neighbours) const {

  int i;

  for (i=0;i<(int)cell.size();i++) {
    if (neighbours[c.sigma][i]==EMPTY) 
      break;
    else
      if (neighbours[c.sigma][i]>0)
	return false;
  }
  
  return true;
  
}



// Based on code by Paulien Hogeweg.
void Dish::CellGrowthAndDivision(void) {
  
  vector<bool> which_cells(cell.size());

  static int mem_area=0;
  
  // if called for the first time: calculate mem_area
  if (!mem_area) {
    mem_area=TargetArea()/CountCells();
  }
  
  int cell_division=0;

  vector<Cell>::iterator c;
  for ( (c=cell.begin(), c++);
	c!=cell.end();
	c++) {
    
    if ( (c->Area()-c->TargetArea())>c->GrowthThreshold() ) {
      c->IncrementTargetArea();
      
    }
    
    if ( (c->TargetArea() > 2 * mem_area ) ) {
      which_cells[c->Sigma()]=true;
      cell_division++;
    }
  }
 
  // Divide scheduled cells
  if (cell_division) {
    CPM->DivideCells(which_cells);
  }

}


int Dish::CountCells(void) const {
  
  int amount=0;
  vector<Cell>::const_iterator i;
  for ( (i=cell.begin(),i++); i!=cell.end(); i++) {
    if (i->AliveP()) {
      amount++;
    } else {
      cerr << "Dead cell\n";
    }
  }
  return amount;
}

 

int Dish::Area(void) const {
  
  int total_area=0;
  
  vector<Cell>::const_iterator i;
  for ( (i=cell.begin(),i++);
	i!=cell.end();
	++i) {
    
    total_area+=i->Area();
    
  }
  return total_area;
}

int Dish::TargetArea(void) const {
  
  int total_area=0;
  
  vector<Cell>::const_iterator i;
  for ( (i=cell.begin(),i++);
	i!=cell.end();
	++i) {
    
    if (i->AliveP()) 
      total_area+=i->TargetArea();
    
  }
  return total_area;
}



void Dish::SetCellOwner(Cell &which_cell) {
  which_cell.owner=this;
}



void Dish::ClearGrads(void) {

  vector<Cell>::iterator i;
  for ( (i=cell.begin(), i++); i!=cell.end(); i++) {
    i->ClearGrad();
  }
}


int Dish::ZygoteArea(void) const {
    return CPM->ZygoteArea();
}

int Dish::Time(void) const {
    return CPM->Time();
}


void Dish::MeasureChemConcentrations(void) {
 
  // clear chemical concentrations
  for (vector<Cell>::iterator c=cell.begin();
       c!=cell.end();
       c++) {
    for (int ch=0;ch<par.n_chem;ch++) 
      c->chem[ch]=0.;
  }

  // calculate current ones
  for (int ch=0;ch<par.n_chem;ch++)
    for (int i=0;i<SizeX()*SizeY();i++) {
      
      int cn=CPM->Sigma(0,i);
      if (cn>=0) 
	cell[cn].chem[ch]+=PDEfield->Sigma(ch,0,i);
	
    }

    for (vector<Cell>::iterator c=cell.begin();
       c!=cell.end();
       c++) {
      for (int ch=0;ch<par.n_chem;ch++) 
	c->chem[ch]/=(double)c->Area();
    }

}

/* void Dish::ExportMultiCellDS (const char *fname) {
    
    FILE *fp=fopen(fname,"w");
    // Export CA field
    CPM->ExportMultiCellDS(fp);
    
    // Export PDE field
    PDEfield->ExportMultiCellDS(fp);
    
    fclose(fp);
    
}*/

void Dish::ExportMultiCellDS (const char *fname) {
    
    // Setup for MultiCellDS output
    MultiCellDS *h = new MultiCellDS;

    // Create a list of individual cells for MultiCellDS
    cell::cell_population_individual *cpi = new cell::cell_population_individual;
 
    //mesh::list_of_voxels = new mesh::list_of_voxels;
    mesh::mesh *mesh = new mesh::mesh;
    
    // get pixels from cells
    vector < list < pair< int, int > > > cellpixels(cell.size());
    
    for (int x=1;x<CPM->SizeX()-1;x++) {
        for (int y=1;y<CPM->SizeY()-1;y++) {
            int s;
            if ((s=CPM->Sigma(x,y))) {
                cellpixels[s].push_back(pair<int, int>(x,y));
            }
        }
    }
    
    // For each individual cell ...
    for (vector<Cell>::iterator c=cell.begin();
         c!=cell.end();
         c++) {
    
        if (c->Sigma()==0) continue;
            // Create a new cell for MultiCellDS
        cell::cell *mcds_cell = new cell::cell;
        cout  << "Created an mcds_cell for cell " << c->Sigma() << "\n";
        // Setup for obtaining different aspects of cellular data
        phenotype_dataset::phenotype_dataset *pds = new phenotype_dataset::phenotype_dataset;
        phenotype_common::geometrical_properties *gp = new phenotype_common::geometrical_properties;
        
        // Translate cell radius from CPM to MultiCellDS
        common::units_decimal *radius_dec = new common::units_decimal; // Setup a units_decimal element
        phenotype_common::lengths *lengths_xml = new phenotype_common::lengths; // Setup up the parent container of radius
        
        cout << "Calculating radius... ";
        double radius = (c->Length()*par.dx)*1e+6; // Take the cell length as radius for now
        radius_dec->units("micron"); // Get and the units
        radius_dec->base_value(radius); // Set the radius
        lengths_xml->radius(radius_dec); // Place the radius in lengths
        gp->lengths(lengths_xml); // Place lengths in geometric properties
        cout << radius << "\n" << endl;
        
        
        // Translate cell volume from CPM to MultiCellDS
        cout << "Calculating volume... ";
        common::units_decimal *volume_dec = new common::units_decimal;
        phenotype_common::volumes *volumes_xml = new phenotype_common::volumes;
        double volume = (c->Area()*par.dx*par.dx)*1e+12; // convert to square microns
        volume_dec->units("square micron");
        volume_dec->base_value(volume);
        volumes_xml->total_volume(volume_dec);
        cout << volume << "\n";

        // Places volumes in geometric properties
        gp->volumes(volumes_xml);
        
        // Create a phenotype
        phenotype::phenotype *p; p = new phenotype::phenotype;
        p->geometrical_properties(gp); // Place geometric properties in phenotype
        phenotype_base::phenotype_type pt = phenotype_base::phenotype_type::current; // Select the current phenotype instead of target or expected
        p->type(pt); // Place the phenotype type into the phenotype
        
        // Put the phenotype into the phenotype dataset
        for(int i=0; i<1; i++)
            pds->phenotype().push_back(p);

        // Add the phenotype dataset to the cell
        mcds_cell->phenotype_dataset(pds);
        
        // Work on the cell state
        state::state *state = new state::state;
        common::units_string *str = new common::units_string;
        common::units_double_list *udl = new common::units_double_list;
        
        
        // Translate position from CPM to MultiCellDS
        double cx, cy, cz;
        cz=0.; // 2D
        c->GetCentroid(&cx, &cy);

        
        // convert from pixel to micron
        cx=(cx*par.dx)*1e+6;
        cy=(cy*par.dx)*1e+6;
        
        udl->push_back(cx); udl->push_back(cy); udl->push_back(cz);
        udl->units("micron");
        state->position(udl);
        
        // add voxels to cell
        //cell::population_vector *popvec = new cell::population_vector;
        /*state::voxels *voxels = new state::voxels;
        
        for (list< pair<int, int> >::const_iterator v=cellpixels[c->Sigma()].begin();
             v!=cellpixels[c->Sigma()].end();
             v++) {
            voxels->push_back((unsigned int)v->first);
            voxels->push_back((unsigned int)v->second);
            voxels->push_back(0);
        }*/
        
        // add voxels to cell
        //cell::population_vector *popvec = new cell::population_vector;
        
        /*mesh::int_list_xpath *voxels = new mesh::int_list_xpath;
        
        for (list< pair<int, int> >::const_iterator v=cellpixels[c->Sigma()].begin();
             v!=cellpixels[c->Sigma()].end();
             v++) {
            voxels->push_back((unsigned int)v->first);
            voxels->push_back((unsigned int)v->second);
            voxels->push_back(0);
            voxels->grouping_number(3); // How many numbers are grouped to make an index
            voxels->xpath("/MultiCellDS/cellular_information/mesh/voxels/"); // Where to find the voxels that the index uses. Use XPATH to navigate through the XML tree.
            
        }*/
        
        
     
        //cell::population_vector = new cell::population_vector;
        
         mesh::int_list_xpath *voxels = new mesh::int_list_xpath;
         for (list< pair<int, int> >::const_iterator v=cellpixels[c->Sigma()].begin();
         v!=cellpixels[c->Sigma()].end();
         v++) {
         voxels->push_back((unsigned int)v->first);
         voxels->push_back((unsigned int)v->second);
         voxels->push_back(0);
         voxels->grouping_number(3); // How many numbers are grouped to make an index
         voxels->xpath("/MultiCellDS/cellular_information/mesh/voxels/"); // Where to find the voxels that the index uses. Use XPATH to navigate through the XML tree.
         
         }
        state->voxels(voxels);

        cout << "(" << cx << ")" << endl;
        // Add state to the cell
        mcds_cell->state(state);
        
        // Add an ID number to the cell
        mcds_cell->ID(c->Sigma());
        
        // Add the cell to the cell population list
        cpi->cell().push_back(mcds_cell);
        
     
    }
    
    
    // Allow cell populations to have a population of individual cells.
    cell::cell_populations *cps = new cell::cell_populations;
    cps->cell_population(cpi);
    cout << "So far so good" << endl;
    
    // Allow cellular information to have cell populations
    cell::cellular_information *ci = new cell::cellular_information;
    ci->cell_populations(cps);
    cout << "So far so good" << endl;
    /*mesh->voxels(voxels);
    ci->mesh(mesh);*/
    
    // Allow the root MultiCellDS element to have cellular information
    h->cellular_information(ci);
    
    MCDS_type mcds_type;// = new MCDS_type;
    mcds_type.value(MCDS_type::value_type::snapshot_simulation); // This is a simulation snapshot (vs experiemnt or clinical)
    h->type(mcds_type); // Assign the type over to the MultiCellDS element
    h->version("0.5.0"); // State the MultiCellDS version number
    
    //############################################################
    // Adding metadata. Read it in from a file. Easier than
    // declaring all kinds of objects.
    //############################################################
    
    // Parse the metadata file.
    //
    MultiCellDS_paggr MultiCellDS_p;
    xml_schema::document_pimpl doc_p (MultiCellDS_p.root_parser (),
                                      MultiCellDS_p.root_name ());
    MultiCellDS_p.pre ();
    
    // Adding metadata
    doc_p.parse("tst_metadata.xml");
    MultiCellDS* h_meta = MultiCellDS_p.post ();
    
    // Update the created and last_modified_times
    //auto current_time = std::chrono::system_clock::now();
    std::time_t t = std::time(NULL);
    struct tm lt = *std::localtime(&t);
    
    // Handle timezones
    time_t offset_second(lt.tm_gmtoff);
    tm off_t = *std::gmtime(&offset_second);
    short timezone_hour = (off_t.tm_hour+12)%24-12, timezone_minute = off_t.tm_min;
    
    xml_schema::date_time new_times(lt.tm_year, lt.tm_mon, lt.tm_mday, lt.tm_hour, lt.tm_min, lt.tm_sec, timezone_hour, timezone_minute);
    h_meta->metadata().created(new_times); // Copy time over
    h_meta->metadata().last_modified(new_times); // Copy time over
    h->metadata(&h_meta->metadata()); // Copy metadata over. Could also use _clone()
    
    
    // Setup for printing the MultiCellDS file
    MultiCellDS_saggr MultiCellDS_s, MultiCellDS_s2;
    
    // doc_s is intended for printing to the screen
    // doc_s2 is intended for printing to a file
    
    // XSD/e setup for printing
    xml_schema::document_simpl doc_s (MultiCellDS_s.root_serializer (),
                                      MultiCellDS_s.root_name ());
    //xml_schema::document_simpl doc_s2 (MultiCellDS_s2.root_serializer (),
//                                       MultiCellDS_s2.root_name ());
    
    // Initialization
    MultiCellDS_s.pre (*h);
    //MultiCellDS_s2.pre (*h);
    // Print to the screen
    //doc_s.serialize (cout, xml_schema::document_simpl::pretty_print);
    
    // Pick a file to write
    ofstream ofs(fname);
    // Print to the file
    
    doc_s.serialize (ofs, xml_schema::document_simpl::pretty_print);
    ofs << endl; // Add an extra newline at the end of the file
    
    // Cleanup of the writing
    MultiCellDS_s.post ();
    //MultiCellDS_s2.post ();
}

 

int Dish::SizeX(void) { return CPM->SizeX(); }
int Dish::SizeY(void) { return CPM->SizeY(); }	
