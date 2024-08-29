#include "cryosettings.h"

using namespace BioFVM;
using namespace PhysiCell;

Cryocell_Settings cryocell_settings;
Cryocell_User_Parameters cryocell_user_parameters;

bool cryocell_config_dom_initialized =false;
pugi::xml_document cryocell_config_doc;
pugi::xml_node cryocell_config_root;


bool load_Cryocell_config_file( std::string filename )
{
	std::cout << "Using config file " << filename << " ... " << std::endl ; 
	pugi::xml_parse_result result = cryocell_config_doc.load_file( filename.c_str()  );
	
	if( result.status != pugi::xml_parse_status::status_ok )
	{
		std::cout << "Error loading " << filename << "!" << std::endl; 
		return false;
	}
	
	cryocell_config_root = cryocell_config_doc.child("cryocell_settings");
	cryocell_config_dom_initialized = true; 
	
	cryocell_settings.read_from_pugixml(); 
	
	
	cryocell_user_parameters.read_from_pugixml( cryocell_config_root ); 

	return true; 	
}
