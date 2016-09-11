/*!
	\brief Text file input handler

	ALL IDS MUST START WITH AN ALPHABETIC LETTER, FOLLOWED BY LETTERS OR DIGITS
//read file
//parse by libxml2
//some check
//write data structure
//write to a file if need compile
//compile
*/

#ifndef _INPUT_ODE_BEFORE_COMPILE_MASS_ACTION_H_
#define _INPUT_ODE_BEFORE_COMPILE_MASS_ACTION_H_

#include "Input.h"
#include <fstream>

//#define _CUSTOM_PROPENSITY_FUNCTIONS_FILENAME_ "CustomPropensityFunctions.h"
namespace STOCHKIT
{

 template<typename _populationVectorType,
	typename _stoichiometryType,
	typename _propensitiesFunctorType,
	typename _dependencyGraphType>
 class Input_ODE_before_compile_mass_action :
	public Input<_populationVectorType, _stoichiometryType, _propensitiesFunctorType, _dependencyGraphType>
 {
	std::vector<std::string> rightHandSide;  // f
	std::vector<std::string> jacobian;       // df/dy
        std::vector<std::string> jacobian_rate;  // df/dc
        std::vector<double>      rateList;       // c
	bool ODE_ready;
	double volume, NATIMESVOLUME, RTOL, ATOL;
        unsigned int MXSTEPS;
	static const double AvogadroConstant = 6.02214E23;
        typedef Input<_populationVectorType, _stoichiometryType, _propensitiesFunctorType, _dependencyGraphType> Base;
        typedef typename Base::_populationValueType _populationValueType;
        using Base::ParametersList;
        using Base::simpleCalculator;
        using Base::NumberOfSpecies;
        using Base::SpeciesList;
 
 public:
        Input_ODE_before_compile_mass_action(char *xmlFilename):
		Input<_populationVectorType, _stoichiometryType, _propensitiesFunctorType, _dependencyGraphType>(xmlFilename),
		rightHandSide(this->NumberOfSpecies, "RCONST(0)"),
		jacobian(this->NumberOfSpecies*this->NumberOfSpecies,"RCONST(0)"),
		jacobian_rate(this->NumberOfSpecies*this->NumberOfReactions,"RCONST(0)"),
		rateList(this->NumberOfReactions,0.0),
		ODE_ready(false),
		NATIMESVOLUME(1.0)
	{
	}

 protected:
	// calculate rate based on the value stored in parameterslist
        double rateCalculation(std::string equation)
	{
		std::vector<unsigned int> ParametersAffectRate;
		std::vector<unsigned int>::iterator para_it; // iterator of parameters in link graph
		
		ParametersAffectRate = this->ParametersList.analyzeParameterExpression(equation);
	
		bool calculationStatus = false;

		for( para_it = ParametersAffectRate.begin(); para_it < ParametersAffectRate.end(); ++para_it ){
			if( this->ParametersList[*para_it].CalculateFlag == -1 ){
				calculationStatus = this->ParametersList.calculateParameter(*para_it);
				if(!calculationStatus){                   
					std::cerr << "StochKit ERROR (Input_mass_action::rateCalculation): while calculating rate " << equation << std::endl;
					return BADRESULT;
				}
			}
		}
		
		std::string substitutedEquation = this->ParametersList.parameterSubstitution(equation);
		if( substitutedEquation.empty() ){
			std::cerr << "StochKit ERROR (Input_mass_action::rateCalculation): while calculating rate " << equation << std::endl;
			return BADRESULT;
		}
		
		return this->simpleCalculator.calculateString(substitutedEquation);
	}

	bool getODEReady()
	{
                double rate;

                typename Input<_populationVectorType, _stoichiometryType, _propensitiesFunctorType, _dependencyGraphType>::Reaction *cur_reaction;

		int volume_position;
		if((volume_position = findParameterWithId("Volume")) != -1){
		} else if((volume_position = findParameterWithId("VOLUME")) != -1){
		} else if((volume_position = findParameterWithId("volume")) != -1){
		} else {
			std::cout << "StochKit WARNING (Input_ODE_before_compile_mass_action::getODEReady): no volume parameter is specified, assuming model file is written in concentration style." << std::endl;
			volume = 1.0;
			NATIMESVOLUME = 1.0;
		}

		if(volume_position != -1){
			if(this->ParametersList[volume_position].CalculateFlag != 1){
				std::cerr << "StochKit ERROR (Input_ODE_before_compile_mass_action::getODEReady): ParametersList not calculated." << std::endl;
				exit(1);
			}
			volume = this->ParametersList[volume_position].Value;
			NATIMESVOLUME = volume * AvogadroConstant;
		}

                std::string rateExp;
		std::string customized_propensities;
		std::ostringstream temp_expressions;

                for(int i=0; i<this->NumberOfReactions; ++i){
                        cur_reaction = &this->ReactionsList[i];

                        if(cur_reaction->Type == 0){
                                rate = rateCalculation(cur_reaction->Rate);
                                if( rate == BADRESULT ){
                                        std::cerr << "StochKit ERROR (Input_ODE_before_compile_mass_action::getODEReady): while calculating rate of reaction " << cur_reaction->Id << "\n";
                                        return false;
                                }
				rateList[i] = rate;
				temp_expressions.str("");
				temp_expressions << "data->p[" << i << "]";
				rateExp = temp_expressions.str();
                                switch ( cur_reaction->Reactants.size() ){
                                        case 0:
//						std::cout << "0-th order: " << i << std::endl;
						for(std::size_t j=0; j<cur_reaction->Products.size(); ++j){
							temp_expressions.str("");
							temp_expressions << " +RCONST(" << (double)cur_reaction->Products[j].Stoichiometry/NATIMESVOLUME << ")*" << rateExp;
							rightHandSide[cur_reaction->Products[j].Index] += temp_expressions.str();
							
							temp_expressions.str("");
							temp_expressions << " +RCONST(" << (double)cur_reaction->Products[j].Stoichiometry/NATIMESVOLUME << ")";
							jacobian_rate[cur_reaction->Products[j].Index * this->NumberOfReactions + i] += temp_expressions.str();
						}
                                                break;
                                        case 1:
                                                if( cur_reaction->Reactants[0].Stoichiometry == -1 ) {
//							std::cout << "1-st order: " << i << std::endl;
							temp_expressions.str("");
							temp_expressions << " +RCONST(" << (double)cur_reaction->Reactants[0].Stoichiometry << ")*" << rateExp << "*Ith(y," << cur_reaction->Reactants[0].Index+1 << ")";
							rightHandSide[cur_reaction->Reactants[0].Index] += temp_expressions.str();
							
							temp_expressions.str("");
							temp_expressions << " +RCONST(" << (double)cur_reaction->Reactants[0].Stoichiometry <<")*" << rateExp ;
							jacobian[cur_reaction->Reactants[0].Index * this->NumberOfSpecies + cur_reaction->Reactants[0].Index] += temp_expressions.str();
							
							temp_expressions.str("");
							temp_expressions << " +RCONST(" << (double)cur_reaction->Reactants[0].Stoichiometry << ")*Ith(y," << cur_reaction->Reactants[0].Index+1 << ")";
							jacobian_rate[cur_reaction->Reactants[0].Index * this->NumberOfReactions + i] += temp_expressions.str();

							for(std::size_t j=0; j<cur_reaction->Products.size(); ++j){
								temp_expressions.str("");
								temp_expressions << " +RCONST(" << (double)cur_reaction->Products[j].Stoichiometry << ")*" << rateExp << "*Ith(y," << cur_reaction->Reactants[0].Index+1 << ")";
								rightHandSide[cur_reaction->Products[j].Index] += temp_expressions.str();
							
								temp_expressions.str("");
								temp_expressions << " +RCONST(" << (double)cur_reaction->Products[j].Stoichiometry <<")*" << rateExp;
								jacobian[cur_reaction->Products[j].Index * this->NumberOfSpecies + cur_reaction->Reactants[0].Index] += temp_expressions.str();
							
								temp_expressions.str("");
								temp_expressions << " +RCONST(" << (double)cur_reaction->Products[j].Stoichiometry << ")*Ith(y," << cur_reaction->Reactants[0].Index+1 << ")";
								jacobian_rate[cur_reaction->Products[j].Index * this->NumberOfReactions + i] += temp_expressions.str();
							}
						}
                                                else if( cur_reaction->Reactants[0].Stoichiometry == -2 ) {
//							std::cout << "2-nd order of same species: " << i << std::endl;
							temp_expressions.str("");
							temp_expressions << " +RCONST(" << (double)cur_reaction->Reactants[0].Stoichiometry*NATIMESVOLUME << ")*" << rateExp << "*Ith(y," << cur_reaction->Reactants[0].Index+1 << ")*Ith(y," << cur_reaction->Reactants[0].Index+1 << ")";
							rightHandSide[cur_reaction->Reactants[0].Index] += temp_expressions.str();
							
							temp_expressions.str("");
							temp_expressions << " +RCONST(" << (double)cur_reaction->Reactants[0].Stoichiometry*2*NATIMESVOLUME << ")*" << rateExp << "*Ith(y," << cur_reaction->Reactants[0].Index+1 << ")";
							jacobian[cur_reaction->Reactants[0].Index * this->NumberOfSpecies + cur_reaction->Reactants[0].Index] += temp_expressions.str();
							
							temp_expressions.str("");
							temp_expressions << " +RCONST(" << (double)cur_reaction->Reactants[0].Stoichiometry*NATIMESVOLUME << ")*Ith(y," << cur_reaction->Reactants[0].Index+1 << ")*Ith(y," << cur_reaction->Reactants[0].Index+1 << ")";
							jacobian_rate[cur_reaction->Reactants[0].Index * this->NumberOfReactions + i] += temp_expressions.str();

							for(std::size_t j=0; j<cur_reaction->Products.size(); ++j){
								temp_expressions.str("");
								temp_expressions << " +RCONST(" << (double)cur_reaction->Products[j].Stoichiometry*NATIMESVOLUME << ")*" << rateExp << "*Ith(y," << cur_reaction->Reactants[0].Index+1 << ")*Ith(y," << cur_reaction->Reactants[0].Index+1 << ")";
								rightHandSide[cur_reaction->Products[j].Index] += temp_expressions.str();
							
								temp_expressions.str("");
								temp_expressions << " +RCONST(" << (double)cur_reaction->Products[j].Stoichiometry*2*NATIMESVOLUME << ")*" << rateExp << "*Ith(y," << cur_reaction->Reactants[0].Index+1 << ")";
								jacobian[cur_reaction->Products[j].Index * this->NumberOfSpecies + cur_reaction->Reactants[0].Index] += temp_expressions.str();
							
								temp_expressions.str("");
								temp_expressions << " +RCONST(" << (double)cur_reaction->Products[j].Stoichiometry*NATIMESVOLUME << ")*Ith(y," << cur_reaction->Reactants[0].Index+1 << ")*Ith(y," << cur_reaction->Reactants[0].Index+1 << ")";
								jacobian_rate[cur_reaction->Products[j].Index * this->NumberOfReactions + i] += temp_expressions.str();
							}
						}
                                                else{
                                                        std::cerr << "StochKit ERROR (Input_ODE_before_compile_mass_action::getODEReady): currently the highest order mass-action reaction supported is bi-molecular reaction while Reaction " << cur_reaction->Id << " is not\n";
                                                        return false;
                                                }
                                                break;
                                        case 2:
                                                if( cur_reaction->Reactants[0].Stoichiometry != -1 || cur_reaction->Reactants[1].Stoichiometry != -1 ){
                                                        std::cerr << "StochKit ERROR (Input_ODE_before_compile_mass_action::getODEReady): currently the highest order mass-action reaction supported is bi-molecular reaction while Reaction " << cur_reaction->Id << " is not\n";
                                                        return false;
						}
						else{
//							std::cout << "2-nd order of different species: " << i << std::endl;
							for(std::size_t j=0; j<cur_reaction->Reactants.size(); ++j){
								temp_expressions.str("");
								temp_expressions << " +RCONST(" << (double)cur_reaction->Reactants[j].Stoichiometry*NATIMESVOLUME << ")*" << rateExp << "*Ith(y," << cur_reaction->Reactants[0].Index+1 << ")*Ith(y," << cur_reaction->Reactants[1].Index+1 << ")";
								rightHandSide[cur_reaction->Reactants[j].Index] += temp_expressions.str();
							
								temp_expressions.str("");
								temp_expressions << " +RCONST(" << (double)cur_reaction->Reactants[j].Stoichiometry*NATIMESVOLUME << ")*" << rateExp << "*Ith(y," << cur_reaction->Reactants[1].Index+1 << ")";
								jacobian[cur_reaction->Reactants[j].Index * this->NumberOfSpecies + cur_reaction->Reactants[0].Index] += temp_expressions.str();
							
								temp_expressions.str("");
								temp_expressions << " +RCONST(" << (double)cur_reaction->Reactants[j].Stoichiometry*NATIMESVOLUME<< ")*" << rateExp << "*Ith(y," << cur_reaction->Reactants[0].Index+1 << ")";
								jacobian[cur_reaction->Reactants[j].Index * this->NumberOfSpecies + cur_reaction->Reactants[1].Index] += temp_expressions.str();
							
								temp_expressions.str("");
								temp_expressions << " +RCONST(" << (double)cur_reaction->Reactants[j].Stoichiometry*NATIMESVOLUME << ")*Ith(y," << cur_reaction->Reactants[0].Index+1 << ")*Ith(y," << cur_reaction->Reactants[1].Index+1 << ")";
								jacobian_rate[cur_reaction->Reactants[j].Index * this->NumberOfReactions + i] += temp_expressions.str();
							}
							for(std::size_t j=0; j<cur_reaction->Products.size(); ++j){
								temp_expressions.str("");
								temp_expressions << " +RCONST(" << (double)cur_reaction->Products[j].Stoichiometry*NATIMESVOLUME<< ")*" << rateExp << "*Ith(y," << cur_reaction->Reactants[0].Index+1 << ")*Ith(y," << cur_reaction->Reactants[1].Index+1 << ")";
								rightHandSide[cur_reaction->Products[j].Index] += temp_expressions.str();
							
								temp_expressions.str("");
								temp_expressions << " +RCONST(" << (double)cur_reaction->Products[j].Stoichiometry*NATIMESVOLUME<< ")*" << rateExp << "*Ith(y," << cur_reaction->Reactants[1].Index+1 << ")";
								jacobian[cur_reaction->Products[j].Index * this->NumberOfSpecies + cur_reaction->Reactants[0].Index] += temp_expressions.str();
							
								temp_expressions.str("");
								temp_expressions << " +RCONST(" << (double)cur_reaction->Products[j].Stoichiometry*NATIMESVOLUME<< ")*" << rateExp << "*Ith(y," << cur_reaction->Reactants[0].Index+1 << ")";
								jacobian[cur_reaction->Products[j].Index * this->NumberOfSpecies + cur_reaction->Reactants[1].Index] += temp_expressions.str();
							
								temp_expressions.str("");
								temp_expressions << " +RCONST(" << (double)cur_reaction->Products[j].Stoichiometry*NATIMESVOLUME<< ")*Ith(y," << cur_reaction->Reactants[0].Index+1 << ")*Ith(y," << cur_reaction->Reactants[1].Index+1 << ")";
								jacobian_rate[cur_reaction->Products[j].Index * this->NumberOfReactions + i] += temp_expressions.str();
							}
                                                }
                                                break;
                                        default:
                                                std::cerr << "StochKit ERROR (Input_ODE_before_compile_mass_action::getODEReady): more than 2 reactants in mass-action reaction " << cur_reaction->Id << "\n";
                                                return false;
                                }
                        } else{
                                  std::cerr << "StochKit ERROR (Input_ODE_before_compile_mass_action::getODEReady): reaction " << cur_reaction->Id << " is not a mass-action\n";
                                  return false;
                        }
                }

		ODE_ready = true;
		return true;
	}

	// substitute parameters and species in custom propensities
	std::string customPropensitySubstitution(std::string equation)
	{
		// locate a parameter name in equation
		unsigned int begin = 0, end = 0;
		std::string substitutedEquation = equation;

		while( begin < substitutedEquation.length() ){
			if( (isalpha(substitutedEquation.at(begin)) || substitutedEquation.at(begin) == '_') && ( (begin==0) || ((substitutedEquation.at(begin) != 'e') && (substitutedEquation.at(begin) != 'E')) || !(isalnum(substitutedEquation.at(begin-1)) || substitutedEquation.at(begin-1) == '_') )){
				end = begin+1;
				while( (end<substitutedEquation.length()) && (isalnum(substitutedEquation.at(end)) || substitutedEquation.at(end) == '_') )
					++end;
				
				std::string parameterName = substitutedEquation.substr(begin,end-begin);
				
				// search for the parameter in ParametersList
				unsigned int i = 0;
				while( (i < this->ParametersList.size()) && (parameterName.compare(this->ParametersList[i].Id)!=0) )
					++i;

				if( i != this->ParametersList.size() ){
					// substitute parameter with its value
//					std::ostringstream parameterValue;
//					parameterValue <<  this->ParametersList[i].Value;
//					substitutedEquation.replace(begin, end-begin, parameterValue.str());
					std::ostringstream parameterValue;
					parameterValue <<  "data->p[" << i << "]";
					substitutedEquation.replace(begin, end-begin, parameterValue.str());
					
					begin = begin + parameterValue.str().size();
				}
				else{
					// search for species in SpeciessList
					unsigned int j = 0;
					while( (j < this->SpeciesList.size()) && (parameterName.compare(this->SpeciesList[j].Id)!=0) )
						++j;

					if( j != this->SpeciesList.size() ){
						// substitute species name with its population variable
						std::ostringstream speciesReference;
						speciesReference << "Ith(y," << j+1 <<")*"<<NATIMESVOLUME;
						substitutedEquation.replace(begin, end-begin, speciesReference.str());
						
						begin = begin + speciesReference.str().size();
					}
					else{
						std::cout << "StochKit WARNING (Input_ODE_before_compile_mass_action::customPropensitySubstitution): function \"" << parameterName << "\" written into custom propensity function, please make sure it's a legitimate c++ function \n";
						begin = begin + parameterName.size();
					}
				}
			}
			else if (substitutedEquation.at(begin) == '/'){
				substitutedEquation.replace(begin, 1, "/(double)");
				begin += 9;
			}
			else{
				++begin;
			}
		}

		return substitutedEquation;
	}

 public:
	bool writeODEFile(char *ODETemplateFileName, char *ODEFileName, std::vector<std::size_t> output_species_index, std::vector<std::string> output_species_names, std::vector<std::string> rate_constants_subset, double RTOL, double ATOL, unsigned int MSTEPS)
	{
		if(!ODE_ready){
			getODEReady();
		}

		std::size_t num_of_par;
		std::vector<std::size_t> rate_constants_index;

		if(rate_constants_subset.size() == 0){
			num_of_par = this->NumberOfReactions;
			rate_constants_index.resize(num_of_par);
	                for (std::size_t i=0; i<num_of_par; ++i) {
				rate_constants_index[i] = i;
			}
		} else{
			num_of_par = rate_constants_subset.size();
			rate_constants_index.resize(num_of_par);
	                for (std::size_t i=0; i<num_of_par; ++i) {
				rate_constants_index[i] = -1;
	                        std::istringstream iss(rate_constants_subset[i]);
	                        std::size_t index;
	                        iss >> index;
	                        if (!iss.fail()) {
                                        std::cout << "StochKit ERROR (Input_ODE_before_compile_mass_action::writeODEFile): rate_constant \""<<rate_constants_subset[i]<<"\" should not be an integer, full reaction name only.\n";
                                        exit(1);
                                } else {
                                	for (std::size_t j=0; j<this->NumberOfReactions; ++j) {
                                        	if (rate_constants_subset[i].compare(this->ReactionsList[j].Id)==0) {
                                               		rate_constants_index[i] = j;
                                                	break;
                                        	}
                                	}
					if(rate_constants_index[i] == -1){
                                        	std::cout << "StochKit ERROR (Input_ODE_before_compile_mass_action::writeODEFile): rate_constant \""<<rate_constants_subset[i]<<"\" can not be found in reactions' list.\n";
                                        	exit(1);
					}
        	                }
	                }
		}
		
		if(output_species_index.size() != output_species_names.size()){
			std::cerr << "StochKit ERROR (Input_ODE_before_compile_mass_action::writeODEFile): output_species_index.size() != output_species_names.size(), something is seriously wrong." << std::endl;
			exit(1);
		} else if(output_species_index.size() == 0){
			std::cerr << "StochKit ERROR (Input_ODE_before_compile_mass_action::writeODEFile): output_species_index.size() == 0. Please set approapriate species to output." << std::endl;
			exit(1);
		}

		std::ifstream ODETemplateFile;
		std::ofstream ODEFile;
		ODETemplateFile.open(ODETemplateFileName);
		ODEFile.open(ODEFileName);

		if( !ODETemplateFile.is_open() ){
			std::cerr << "StochKit ERROR (Input_ODE_before_compile_mass_action::writeODEFile): open ODE template file " << ODETemplateFileName << " failed." << std::endl;
			exit(1);
		} else if( !ODEFile.is_open() ) {
			std::cerr << "StochKit ERROR (Input_ODE_before_compile_mass_action::writeODEFile): open write-to ODE file " << ODEFileName << " failed." << std::endl;
			exit(1);
		}

		_populationVectorType initialPop = writeInitialConcentration();

		std::string line;
		std::size_t insert_found;

		while(ODETemplateFile.good()){
			getline(ODETemplateFile, line);
			insert_found = line.find("INSERT GENERATED");
			if(insert_found == std::string::npos){
				ODEFile << line << std::endl;
			} else {
				if(line.find("CONSTANTS") != std::string::npos){
					ODEFile << line << std::endl;
					ODEFile << "#define RTOL  RCONST(" << RTOL << ")   /* scalar relative tolerance            */" << std::endl;
					ODEFile << "#define ATOL  RCONST(" << ATOL << ")    /* scalar absolute tolerance components */" << std::endl;
					ODEFile << "#define MXSTEPS  " << MXSTEPS << "           /* max steps before tout */" << std::endl;
					ODEFile << "#define NEQ   " << this->NumberOfSpecies << "          /* number of equations  */" << std::endl;
					ODEFile << "#define OUTPUT_SPECIES_NUMBER  " << output_species_index.size() << "    /* number of species in output  */" << std::endl;
					ODEFile << "#define NUMPAR  " << this->NumberOfReactions << "    /* number of parameters/rate constants */" << std::endl;
					ODEFile << "#define NUMPAR_SENSI  " << num_of_par << "    /* number of parameters/rate constants in sensitivity analysis */" << std::endl;
				} else if (line.find("INITIAL CONDITION") != std::string::npos){
					ODEFile << line << std::endl;
					for(int i=0;i<this->NumberOfSpecies;++i){
						ODEFile << "  Ith(y," << i+1 << ") = RCONST(" << (double)initialPop[i]/volume <<");" << std::endl;
					}
				} else if (line.find("PARAMETERS LIST") != std::string::npos){
					ODEFile << line << std::endl;
					for(int i=0;i<this->NumberOfReactions;++i){
						ODEFile << "  data->p[" << i << "] = RCONST(" << rateList[i] <<");" << std::endl;
					}
				} else if (line.find("RIGHT HAND SIDE") != std::string::npos){
					ODEFile << line << std::endl;
					for(int i=0;i<this->NumberOfSpecies;++i){
						ODEFile << "  Ith(ydot," << i+1 << ") = " << rightHandSide[i] << ";" << std::endl;
					}
				} else if (line.find("USER SUPPLIED JACOBIAN") != std::string::npos){
					ODEFile << line << std::endl;
					for(int i=0;i<this->NumberOfSpecies;++i){
						for(int j=0;j<this->NumberOfSpecies;++j){
							ODEFile << "  IJth(J," << i+1 << "," << j+1 <<") = " << jacobian[i*this->NumberOfSpecies+j] << ";" << std::endl;
						}
					}
				} else if (line.find("USER SUPPLIED FORWARD SENSITIVITY") != std::string::npos){
					ODEFile << line << std::endl;
					for(int i=0;i<this->NumberOfSpecies;++i){
						ODEFile << "  sd[" << i << "] = ";
                                                for(int j=0;j<this->NumberOfSpecies;++j){
							ODEFile << "+(" << jacobian[i*this->NumberOfSpecies+j] << ")*s[" << j << "]";
						}
						ODEFile << ";" << std::endl;
					}
					ODEFile << "  switch (iS) {" << std::endl;
					for(int j=0;j<num_of_par;++j){
						std::size_t index = rate_constants_index[j];
						ODEFile << "   case "<< j << ":" << std::endl;
					        for(int i=0;i<this->NumberOfSpecies;++i){
							ODEFile << "      sd[" << i << "] += " << jacobian_rate[i * this->NumberOfReactions + index] << ";" << std::endl;
						}
						ODEFile << "   break;" << std::endl;
					}
					ODEFile << "  }" << std::endl;
				} else if (line.find("OUTPUT SPECIES NAMES") != std::string::npos){
					ODEFile << line << std::endl;
					ODEFile << "  char *species_names[" << output_species_index.size() << "];" << std::endl;
					for(std::size_t i=0;i<output_species_index.size();++i){
						ODEFile << "   species_names[" << i << "] = \"" << output_species_names[i] << "\";" << std::endl;
					}
				} else if (line.find("OUTPUT SPECIES INDEX") != std::string::npos){
					ODEFile << line << std::endl;
					ODEFile << "  int species_index[" << output_species_index.size() << "] = { ";
					for(std::size_t i=0;i<output_species_index.size()-1;++i){
						ODEFile << " " << output_species_index[i] << ", ";
					}
					ODEFile << " " << output_species_index[output_species_index.size()-1] << " ";
					ODEFile << " };" << std::endl;
				} else if (line.find("PLIST") != std::string::npos){
					ODEFile << line << std::endl;
					ODEFile << "  int plist[NUMPAR_SENSI] = {";
					for(int j=0;j<num_of_par-1;++j){
						ODEFile << " " << rate_constants_index[j] << ",";
					}
					ODEFile << " " << rate_constants_index[num_of_par-1];
					ODEFile << "};" << std::endl;
				} else if (line.find("OUTPUT PARAMETERS NAMES") != std::string::npos){
					ODEFile << line << std::endl;
					ODEFile << "  char *parameters_names[NUMPAR_SENSI];" << std::endl;
					for(int j=0;j<num_of_par;++j){
						ODEFile << "   parameters_names[" << j << "] = \"" << this->ReactionsList[rate_constants_index[j]].Id << "\";" << std::endl;
					}
				} else{
					std::cerr << "StochKit ERROR (Input_ODE_before_compile_mass_action::writeODEFile): template file corrupted, indication line not recognized: " << line;
					ODETemplateFile.close();
					ODEFile.close();
					exit(1);
				}
/*
				while(ODETemplateFile.good()){
					getline(ODETemplateFile, line);
					end_found = line.find("END OF GENERATED");
					if(end_found != std::string::npos){
						ODEFile << line;
						break;
					}
				}
*/
			}
		}

		ODETemplateFile.close();
		ODEFile.close();
		
		return true;
	}

  protected:
	int findParameterWithId(std::string givenId)
	{
	        int found = -1;
	        for( unsigned int i = 0; i < this->ParametersList.size(); ++i ){
        	        if(this->ParametersList[i].Id.compare(givenId) == 0){
                	        found = (int)i;
                        	break;
	                }
        	}
	        return found;
	}

_populationValueType concentrationCalculation(std::string equation)
{
        std::vector<unsigned int> ParametersAffectRate;
        std::vector<unsigned int>::iterator para_it; // iterator of parameters in link graph

        ParametersAffectRate = ParametersList.analyzeParameterExpression(equation);

        bool calculationStatus = false;

        for( para_it = ParametersAffectRate.begin(); para_it < ParametersAffectRate.end(); ++para_it ){
                if( ParametersList[*para_it].CalculateFlag == -1 ){
                        calculationStatus = ParametersList.calculateParameter(*para_it);
                        if(!calculationStatus){
                                std::cerr << "StochKit ERROR (Input::populationCalculation): while calculating rate " << equation << std::endl;
                                return BADRESULT;
                        }
                }
        }

        std::string substitutedEquation = customPropensitySubstitution(equation);
        if( substitutedEquation.empty() ){
                std::cerr << "StochKit ERROR (Input::populationCalculation): while calculating initial population " << equation << std::endl;
                return BADRESULT;
        }

        double calculatedPopulation = simpleCalculator.calculateString(substitutedEquation);

        return (_populationValueType)calculatedPopulation;
}

_populationVectorType writeInitialConcentration()
{
          _populationValueType cur_population;
          _populationVectorType X(NumberOfSpecies);

          for(int i=0; i<NumberOfSpecies; ++i){
                cur_population = concentrationCalculation(SpeciesList[i].InitialPopulation);
                if( cur_population == BADRESULT ){
                        std::cerr << "StochKit ERROR (Input_ODE_before_compilation::writeInitialConcentration): while calculating initial population of " << SpeciesList[i].Id << std::endl;
                        exit(1);
                }

                X[i] = cur_population;
          }

          return X;
}

 };

}
#endif

