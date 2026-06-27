#pragma once

#include <algorithm>
#include <cctype>
#include <fstream>
#include <set>
#include <sstream>
#include <string>

#include <nlohmann/json.hpp>

#include <TNL/Config/ConfigDescription.h>
#include <TNL/Config/ConfigEntry.h>
#include <TNL/Config/ConfigEntryList.h>
#include <TNL/Config/ParameterContainer.h>
#include <TNL/Config/parseCommandLine.h>  // for addDefaultValues and checkEnumValues
#include <TNL/Config/parseINIConfigFile.h>  // for the INI fallback in the dispatcher

namespace TNL::SPH {

/**
 * \brief Parses a JSONC configuration file (JSON with comments) into a
 *        \ref TNL::Config::ParameterContainer.
 *
 * This function is a drop-in replacement for \ref TNL::Config::parseINIConfigFile.
 * It reads a flat JSON object whose keys match the entries defined in the
 * \ref TNL::Config::ConfigDescription, converts each JSON value to the C++ type
 * declared in the schema, and populates the returned container.  Both line
 * comments (double-slash) and block comments (slash-asterisk) are supported
 * through nlohmann/json's \c ignore_comments option.
 *
 * \param configPath  Path to the JSONC file.
 * \param description Schema describing expected entries, their types, default
 *                    values and enum constraints.
 * \returns A populated \ref TNL::Config::ParameterContainer.
 * \throws std::runtime_error if the file cannot be opened, the top-level JSON
 *         value is not an object, or an unknown option is encountered.
 * \throws nlohmann::json::exception on JSON parse or type-conversion errors.
 */
[[nodiscard]] inline TNL::Config::ParameterContainer
parseJSONConfigFile( const std::string& configPath, const TNL::Config::ConfigDescription& description )
{
   TNL::Config::ParameterContainer parameters;

   // --- read the file -------------------------------------------------------
   std::ifstream file( configPath );
   if( ! file.is_open() )
      throw std::runtime_error( "Failed to open the configuration file: " + configPath );

   std::stringstream buffer;
   buffer << file.rdbuf();

   // Parse JSONC: the last argument (ignore_comments) enables line and block comments
   nlohmann::json data = nlohmann::json::parse( buffer.str(), nullptr, /*allow_exceptions*/ true, /*ignore_comments*/ true );

   if( ! data.is_object() )
      throw std::runtime_error( "The configuration file " + configPath + " must contain a JSON object at the top level." );

   // --- iterate over all entries in the JSON object -------------------------
   std::set< std::string > undefinedOptions;

   for( auto it = data.begin(); it != data.end(); ++it ) {
      const std::string name = it.key();
      const auto* entryBase = description.getEntry( name );
      if( entryBase == nullptr ) {
         undefinedOptions.insert( name );
         continue;
      }
      const std::string entryType = entryBase->getUIEntryType();
      const auto& value = it.value();

      if( entryType == "bool" ) {
         parameters.addParameter< bool >( name, value.get< bool >() );
         continue;
      }
      else if( entryType == "integer" ) {
         const auto v = value.get< TNL::Config::Integer >();
         const auto& entry = dynamic_cast< const TNL::Config::ConfigEntry< TNL::Config::Integer >& >( *entryBase );
         TNL::Config::checkEnumValues( entry, name, v );
         parameters.addParameter< TNL::Config::Integer >( name, v );
      }
      else if( entryType == "unsigned integer" ) {
         const auto v = value.get< TNL::Config::UnsignedInteger >();
         const auto& entry = dynamic_cast< const TNL::Config::ConfigEntry< TNL::Config::UnsignedInteger >& >( *entryBase );
         TNL::Config::checkEnumValues( entry, name, v );
         parameters.addParameter< TNL::Config::UnsignedInteger >( name, v );
      }
      else if( entryType == "real" ) {
         const auto v = value.get< double >();
         const auto& entry = dynamic_cast< const TNL::Config::ConfigEntry< double >& >( *entryBase );
         TNL::Config::checkEnumValues( entry, name, v );
         parameters.addParameter< double >( name, v );
      }
      else if( entryType == "string" ) {
         const auto v = value.get< std::string >();
         const auto& entry = dynamic_cast< const TNL::Config::ConfigEntry< std::string >& >( *entryBase );
         TNL::Config::checkEnumValues( entry, name, v );
         parameters.addParameter< std::string >( name, v );
         continue;
      }
      else if( entryType == "list of bool" ) {
         const auto list = value.get< std::vector< bool > >();
         parameters.addList< bool >( name, list );
      }
      else if( entryType == "list of integer" ) {
         const auto list = value.get< std::vector< TNL::Config::Integer > >();
         const auto& entry = dynamic_cast< const TNL::Config::ConfigEntry< TNL::Config::Integer >& >( *entryBase );
         TNL::Config::checkEnumValues( entry, name, list );
         parameters.addList< TNL::Config::Integer >( name, list );
      }
      else if( entryType == "list of unsigned integer" ) {
         const auto list = value.get< std::vector< TNL::Config::UnsignedInteger > >();
         const auto& entry = dynamic_cast< const TNL::Config::ConfigEntry< TNL::Config::UnsignedInteger >& >( *entryBase );
         TNL::Config::checkEnumValues( entry, name, list );
         parameters.addList< TNL::Config::UnsignedInteger >( name, list );
      }
      else if( entryType == "list of real" ) {
         const auto list = value.get< std::vector< double > >();
         const auto& entry = dynamic_cast< const TNL::Config::ConfigEntry< double >& >( *entryBase );
         TNL::Config::checkEnumValues( entry, name, list );
         parameters.addList< double >( name, list );
      }
      else if( entryType == "list of string" ) {
         const auto list = value.get< std::vector< std::string > >();
         const auto& entry = dynamic_cast< const TNL::Config::ConfigEntry< std::string >& >( *entryBase );
         TNL::Config::checkEnumValues( entry, name, list );
         parameters.addList< std::string >( name, list );
      }
      else
         // this will not happen if all entry types are handled above
         throw std::runtime_error( "Function parseJSONConfigFile encountered unsupported entry type: " + entryType );
   }

   if( ! undefinedOptions.empty() ) {
      std::string msg = "The configuration file contains the following options which are not defined in the program:\n";
      for( const auto& option : undefinedOptions )
         msg += " - " + option + "\n";
      throw std::runtime_error( msg );
   }

   // add default values for entries that were not present in the file
   TNL::Config::addDefaultValues( description, parameters );

   return parameters;
}

/**
 * \brief Format-agnostic configuration file parser.
 *
 * Dispatches to \ref parseJSONConfigFile for \c .json / \c .jsonc files and to
 * \ref TNL::Config::parseINIConfigFile for all other extensions (including
 * \c .ini).  This allows mixing INI and JSONC configuration files within the
 * same project — each config file is parsed according to its extension.
 *
 * \param configPath  Path to the configuration file.
 * \param description Schema describing expected entries.
 * \returns A populated \ref TNL::Config::ParameterContainer.
 */
[[nodiscard]] inline TNL::Config::ParameterContainer
parseConfigFile( const std::string& configPath, const TNL::Config::ConfigDescription& description )
{
   const auto dotPos = configPath.find_last_of( '.' );
   if( dotPos != std::string::npos ) {
      std::string ext = configPath.substr( dotPos + 1 );
      std::transform( ext.begin(), ext.end(), ext.begin(), []( unsigned char c ) { return static_cast< char >( std::tolower( c ) ); } );
      if( ext == "json" || ext == "jsonc" )
         return parseJSONConfigFile( configPath, description );
   }
   // default: INI
   return TNL::Config::parseINIConfigFile( configPath, description );
}

}  // namespace TNL::SPH
