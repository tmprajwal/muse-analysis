#include <hough.h>

// Instantiate the geometry
void hough::GenerateDet()
{

  // Import the default geometry file or use the command line geometry
  char buffer[1024];
  geofile = "muse_v4"; // Default geometry file

  // Write out the file name
  snprintf(buffer,1000,"%s/.muse/shared/gdml/%s.gdml",getenv("COOKERHOME"),geofile);
  // if(viscotab)
     //gGeoManager=gEve->GetGeometry(buffer);
  // else
    gGeoManager=TGeoManager::Import(buffer);
  
  gGeoManager->DefaultColors();
  gGeoManager->CloseGeometry();

  // // Wrap the geometry for EVE, making sure to incorporate all needed levels
  // etn = new TEveGeoTopNode(gGeoManager,gGeoManager->GetTopNode(),0,7);
  // gEve->AddElement(etn);

  // std::cout<<"\n\nConfiguring initial visualization settings.\n\nIgnore TGeoManager and TGLViewerbase info/errors...\n\n";

  // // Get the pointers to all the nodes you might want to edit
  findnodes();

};


void hough::printallnodes(TGeoNode *node,int level)
{
     TObjArray *k=node->GetNodes();
     TGeoVolume *mat;
     //mat = node->GetPhysVolume();
     //TGeoPhysicalNode *b;
     mat = node->GetVolume();
     int copynr = node->GetNumber();
     TString name = mat->GetName();
      //No one should ever look at this, or ever use this, or ever acknowledge this exists
      //Terrible hack because I do not know Geant4 but I mangled something together
      //So that we could actually have straws in our Event display
      //Please forgive me.
     //std::cout << node->GetName()<<std::endl;
      if(name.BeginsWith("STTL1_plane60"))
      {
        for(int j = 0; j<10; j++)
        {
          TString straw = "STTL1_plane60_"+std::to_string(j);

          int uppernum = 0;
          if(j<5&&(j%2==0))
            uppernum = 55;
          else if(j<5 && j%2==1)
            uppernum = 54;
          else if(j>=5 &&j%2==0)
            uppernum = 54;
          else if(j>=5 &&j%2==1)
            uppernum = 55;

          for(int i = 0; i<uppernum; i++)
          {
            straw = "STTL1_plane60_"+std::to_string(j)+"straw"+std::to_string(i)+"_log";
            if(name.BeginsWith(straw))
            {
              if(uppernum==54)
                STTL160[j][uppernum] = NULL;
              if(j>4)
                STTL160[j][i] = node;
              else
                STTL160[j][uppernum-i-1] = node;
            }
          }
        }
      }
      else if(name.BeginsWith("STTL1_plane90"))
      {
        for(int j = 0; j<10; j++)
        {
          TString straw = "STTL1_plane90_"+std::to_string(j);
          int uppernum = 0;
          if(j<5&&(j%2==0))
            uppernum = 89;
          else if(j<5 && j%2==1)
            uppernum = 88;
          else if(j>=5 &&j%2==0)
            uppernum = 88;
          else if(j>=5 &&j%2==1)
            uppernum = 89;

          for(int i = 0; i<uppernum; i++)
          {
            straw = "STTL1_plane90_"+std::to_string(j)+"straw"+std::to_string(i)+"_log";
            if(name.BeginsWith(straw))
            {
              if(uppernum==88)
                STTL190[j][uppernum] = NULL;
              if(j>4)
                STTL190[j][i] = node;
              else
                STTL190[j][uppernum-i-1] = node;
            }
          }
        }
      }

      //std::cout << STTL1[10]->GetName() << " " << copynr<< " ";
     //std::cout<<level<<" "<<node->GetName()  << "\n";

      if (k) for (int j=0;j<k->GetEntriesFast();j++)
         printallnodes((TGeoNode*) k->At(j),level+1);
 
}


Long_t hough::findnodes()
{

  int all = 0;  // Counter to keep track of if the whole geometry is found
  //std::cout<<"\\n\\nall = "<<all<<"\\n\\n";

  int ind = 0;
  TObjArray *l=gGeoManager->GetTopVolume()->GetNodes();
  printallnodes(gGeoManager->GetTopNode(),0);
      
    
  // Return the bit count (0 if no failure)
  return (Long_t)all;

}