#ifndef __CRYOSETTINGS_H__
#define __CRYOSETTINGS_H__
#include <iostream>
#include <ctime>

#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include <vector>
#include <random>
#include <chrono>
#include <unordered_map>

#include "../../modules/PhysiCell_pugixml.h"
#include "../../BioFVM/BioFVM.h"

#include "../../core/PhysiCell_constants.h" 
#include "../../core/PhysiCell_utilities.h"
using namespace BioFVM;
using namespace PhysiCell;

extern pugi::xml_node cryocell_config_root;

bool load_Cryocell_config_file( std::string filename);



class Cryocell_Settings {
  private:
  public:
  Cryocell_Settings();
  void read_from_pugixml(void);
};

template <class T> 
class Cryo_Parameter
{
 private:
	template <class Y>
	friend std::ostream& operator<<(std::ostream& os, const Cryo_Parameter<Y>& param); 

 public: 
	std::string name; 
	std::string units; 
	T value; 
	
	Cryo_Parameter();
	Cryo_Parameter( std::string my_name ); 
	
	void operator=( T& rhs ); 
	void operator=( T rhs ); 
	void operator=( Cryo_Parameter& p ); 
};

template <class T>
class Cryo_Parameters
{
 private:
	std::unordered_map<std::string,int> cryocell_name_to_index_map; 
	
	template <class Y>
	friend std::ostream& operator<<( std::ostream& os , const Cryo_Parameters<Y>& params ); 

 public: 
	Cryo_Parameters(); 
 
	std::vector< Cryo_Parameter<T> > parameters; 
	
	void add_cryo_parameter( std::string my_name ); 
	void add_cryo_parameter( std::string my_name , T my_value ); 
//	void add_parameter( std::string my_name , T my_value ); 
	void add_cryo_parameter( std::string my_name , T my_value , std::string my_units ); 
//	void add_parameter( std::string my_name , T my_value , std::string my_units ); 
	
	void add_cryo_parameter( Cryo_Parameter<T> param );
	
	int find_cryo_index( std::string search_name ); 
	
	// these access the values 
	T& operator()( int i );
	T& operator()( std::string str ); 

	// these access the full, raw parameters 
	Cryo_Parameter<T>& operator[]( int i );
	Cryo_Parameter<T>& operator[]( std::string str ); 
	
	int size( void ) const; 
};

class Cryocell_User_Parameters
{
 private:
	friend std::ostream& operator<<( std::ostream& os , const Cryocell_User_Parameters up ); 
 
 public:
	Cryo_Parameters<bool> bools; 
	Cryo_Parameters<int> ints; 
	Cryo_Parameters<double> doubles; 
	Cryo_Parameters<std::string> strings; 
	
	void read_from_pugixml( pugi::xml_node parent_node );
}; 


extern Cryocell_Settings cryocell_settings; 

extern Cryocell_User_Parameters cryocell_user_parameters; 



#endif // !__CRYOSETTINGS_H__
