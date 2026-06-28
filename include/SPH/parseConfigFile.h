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
 * \brief Parses a JSONC string (JSON with comments) into a
 *        \ref TNL::Config::ParameterContainer.
 *
 * Core parsing routine used by both \ref parseJSONConfigFile (file input) and
 * inline subsection parsing.  When a string-typed entry receives a JSON object
 * or array as its value, the value is serialized back to a JSON string via
 * \c nlohmann::json::dump().  This enables nested configuration subsections:
 * the caller stores the serialized string and later re-parses it with a
 * sub-schema (see SimulationMonitor for the "measuretool" subsection).
 *
 * \param jsonString  JSONC content as a string.
 * \param sourceLabel  Human-readable label for error messages (e.g. file path).
 * \param description  Schema describing expected entries, their types, default
 *                     values and enum constraints.
 * \returns A populated \ref TNL::Config::ParameterContainer.
 */
[[nodiscard]] inline TNL::Config::ParameterContainer
parseJSONConfigString( const std::string& jsonString,
                       const std::string& sourceLabel,
                       const TNL::Config::ConfigDescription& description )
{
   TNL::Config::ParameterContainer parameters;

   // The last argument (ignore_comments) enables line and block comments
   nlohmann::json data = nlohmann::json::parse( jsonString, nullptr, /*allow_exceptions*/ true, /*ignore_comments*/ true );

   if( ! data.is_object() )
      throw std::runtime_error( "The configuration source " + sourceLabel + " must contain a JSON object at the top level." );

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
         // Nested JSON objects/arrays are serialized to a string so they can be
         // re-parsed later with a sub-schema (enables inline config subsections).
         std::string v;
         if( value.is_object() || value.is_array() )
            v = value.dump();
         else
            v = value.get< std::string >();
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
         throw std::runtime_error( "Function parseJSONConfigString encountered unsupported entry type: " + entryType );
   }

   if( ! undefinedOptions.empty() ) {
      std::string msg = "The configuration source " + sourceLabel + " contains the following options which are not defined in the program:\n";
      for( const auto& option : undefinedOptions )
         msg += " - " + option + "\n";
      throw std::runtime_error( msg );
   }

   TNL::Config::addDefaultValues( description, parameters );

   return parameters;
}

/**
 * \brief Parses a JSONC configuration file (JSON with comments) into a
 *        \ref TNL::Config::ParameterContainer.
 *
 * Reads the file and delegates to \ref parseJSONConfigString.
 *
 * \param configPath  Path to the JSONC file.
 * \param description Schema describing expected entries.
 * \returns A populated \ref TNL::Config::ParameterContainer.
 */
[[nodiscard]] inline TNL::Config::ParameterContainer
parseJSONConfigFile( const std::string& configPath, const TNL::Config::ConfigDescription& description )
{
   std::ifstream file( configPath );
   if( ! file.is_open() )
      throw std::runtime_error( "Failed to open the configuration file: " + configPath );

   std::stringstream buffer;
   buffer << file.rdbuf();

   return parseJSONConfigString( buffer.str(), configPath, description );
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
