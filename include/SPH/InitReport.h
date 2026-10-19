#pragma once

#include <string>

#include <nlohmann/json_fwd.hpp>

namespace TNL {
namespace SPH {

/**
 * \brief Creates the root document of the solver initialization report.
 *
 * The document is pre-seeded with the report schema version and all top-level
 * sections (empty by default). Individual sections are filled by the free
 * builder functions (e.g. \ref addConfigSection) and by the per-class
 * \e writeInitReport methods of the particle sets, boundary patches,
 * multiresolution buffers and decomposition topology, which gather their own
 * state. The whole document is finally emitted with \ref writeInitReportFile.
 *
 * \returns A JSON object with the keys \e schema_version, \e config,
 *          \e fluids, \e boundaries, \e open_boundary_patches,
 *          \e periodic_patches, \e multiresolution_buffers and
 *          \e decomposition.
 */
[[nodiscard]] inline nlohmann::json
makeInitReportDocument();

/**
 * \brief Fills the \e config section of the initialization report from the
 *        solver parameter container.
 *
 * Emits the domain geometry with the same key spellings as the JSONC
 * configuration file: \e searchRadius, \e domainOrigin and \e domainSize
 * (three-component vectors read via the container's \e getXyz helper, exactly
 * as DecompositionTopology.h does), the model resolution parameters \e dp and
 * \e h when the container carries them, and \e numberOfSubdomains (1 when the
 * entry is absent or unset).
 *
 * \tparam Params     Parameter-container-like type (e.g. TNL::Config::ParameterContainer).
 * \param doc         Report document to update (its \e config key is replaced).
 * \param parameters  Solver runtime parameters.
 */
template< typename Params >
void
addConfigSection( nlohmann::json& doc, const Params& parameters );

/**
 * \brief Writes the initialization report document into a JSON file.
 *
 * Writes \e init_report.json with a two-space indent into \e outputDirectory.
 * Under MPI (initialized, more than one process) the file is per-rank:
 * \e init_report_subdomain\<rank\>.json. The directory itself is never
 * created by this function.
 *
 * \param doc              Report document to write.
 * \param outputDirectory  Existing output directory.
 *
 * \throws std::runtime_error when the report file cannot be opened.
 */
inline void
writeInitReportFile( const nlohmann::json& doc, const std::string& outputDirectory );

}  //namespace SPH
}  //namespace TNL

#include "InitReport.hpp"
