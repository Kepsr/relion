/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/
/***************************************************************************
 *
 * Authors:	J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307 USA
 *
 *	All comments concerning this program package may be sent to the
 *	e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef METADATA_LABEL_H
#define METADATA_LABEL_H

#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "src/funcs.h"

class StaticInitialization;

class EMDL {

    public:

    struct LabelData;  // For storing the name and type of a data label

    // This enum defines what MetaDataLabels this class can manage. 
    // If you need a new one add it here and modify affected methods:
    //
    // - static EMDLabel codifyLabel( std::string strLabel ); EMDL::addLabel<int>(EMDL::OPTIMISER_RANDOM_SEED, "randomSeed");

    // - static std::string EMDL::label2Str( EMDLabel inputLabel );
    // - void writeValuesToFile( std::ofstream &outfile, EMDLabel inputLabel );
    // - void addValue( std::string name, std::string value );
    //
    // Keep this special structure (using EMDL::FIRSTLABEL and EMDL::LAST_LABEL)
    // so the programmer can iterate through it like this:
    //
    //	for (EMDLabel mdl = EMDL::FIRST_LABEL; mdl < EMDL::LAST_LABEL; ++mdl)
    //
    enum EMDLabel {
        UNDEFINED = -1, // Keep the order the same as in StaticInitialization below!!
        FIRST_LABEL,
        OBJID = FIRST_LABEL, ///< object id (int), NOTE: This label is special and shouldn't be used

        AREA_ID,   ///< ID   for the area (or field of view). If one does not use (tilt) series, area would be the same as micrograph...
        AREA_NAME, ///< Name for the area (or field of view). If one does not use (tilt) series, area would be the same as micrograph...
        COMMENT, // The COMMENT is handled specially as well

        BODY_MASK_NAME,  ///< For multi-body refinements
        BODY_KEEP_FIXED, ///< For multi-body refinements
        BODY_REFERENCE_NAME,
        BODY_ROTATE_DIRECTION_X,
        BODY_ROTATE_DIRECTION_Y,
        BODY_ROTATE_DIRECTION_Z,
        BODY_ROTATE_RELATIVE_TO,
        BODY_SIGMA_ANG,
        BODY_SIGMA_OFFSET, // deprecated
        BODY_SIGMA_OFFSET_ANGSTROM,
        BODY_SIGMA_ROT,
        BODY_SIGMA_TILT,
        BODY_SIGMA_PSI,
        BODY_STAR_FILE,

        CTF_ASTIGMATISM,
        CTF_BFACTOR,                   ///< B-factor
        CTF_MAXRES,                    ///< Maximum resolution with Thon rings
        CTF_VALIDATIONSCORE,           ///< Gctf-based validation score for CTF fit
        CTF_SCALEFACTOR,               ///< linear scale-factor
        CTF_SAMPLING_RATE,             ///< Sampling rate
        CTF_VOLTAGE,                   ///< Microscope voltage (kV)
        CTF_DEFOCUSU,                  ///< Defocus U (Angstroms)
        CTF_DEFOCUSV,                  ///< Defocus V (Angstroms)
        CTF_DEFOCUS_ANGLE,             ///< Defocus angle (degrees)
        CTF_CS,                        ///< Spherical aberration
        CTF_CA,                        ///< Chromatic aberration
        CTF_DETECTOR_PIXEL_SIZE,       ///< Pixel size for detector as used in CTF-determination (deprecated)
        CTF_POWER_SPECTRUM,
        CTF_ENERGY_LOSS,               ///< Energy loss
        CTF_FOM,                       ///< ctffind FOM (CC) for quality of CTF-fit
        CTF_IMAGE,                     ///< name of an image describing the CTF model
        CTF_LENS_STABILITY,            ///< Lens stability
        CTF_MAGNIFICATION,             ///< Magnification used for CTF-determination (deprecated)
        CTF_PHASESHIFT,                ///< Phase-shift from a phase plate
        CTF_CONVERGENCE_CONE,          ///< Convergence cone
        CTF_LONGITUDINAL_DISPLACEMENT, ///< Longitudinal displacement
        CTF_TRANSVERSAL_DISPLACEMENT,  ///< Transverse   displacement
        CTF_Q0,                        ///< Amplitude contrast
        CTF_K,                         ///< CTF gain
        CTF_VALUE,                     ///< CTF value

        IMAGE_NAME,
        IMAGE_ORI_NAME,
        IMAGE_RECONSTRUCT_NAME,
        IMAGE_ID,
        IMAGE_ENABLED,
        IMAGE_DATATYPE,
        IMAGE_DIMENSIONALITY,
        IMAGE_BEAMTILT_X,
        IMAGE_BEAMTILT_Y,
        IMAGE_MTF_FILENAME,
        IMAGE_OPTICS_GROUP,
        IMAGE_OPTICS_GROUP_NAME,
        IMAGE_ODD_ZERNIKE_COEFFS,
        IMAGE_EVEN_ZERNIKE_COEFFS,
        IMAGE_PIXEL_SIZE,
        IMAGE_MAG_MATRIX_00,
        IMAGE_MAG_MATRIX_01,
        IMAGE_MAG_MATRIX_10,
        IMAGE_MAG_MATRIX_11,

        IMAGE_COORD_X,
        IMAGE_COORD_Y,
        IMAGE_COORD_Z,
        IMAGE_FRAME_NR,
        IMAGE_MAGNIFICATION_CORRECTION,
        IMAGE_NORM_CORRECTION,
        IMAGE_SAMPLINGRATE,
        IMAGE_SAMPLINGRATE_X,
        IMAGE_SAMPLINGRATE_Y,
        IMAGE_SAMPLINGRATE_Z,
        IMAGE_SIZE,
        IMAGE_SIZE_X,
        IMAGE_SIZE_Y,
        IMAGE_SIZE_Z,
        IMAGE_STATS_MIN,
        IMAGE_STATS_MAX,
        IMAGE_STATS_AVG,
        IMAGE_STATS_STDDEV,
        IMAGE_STATS_SKEW,
        IMAGE_STATS_KURT,
        IMAGE_WEIGHT,

        JOB_IS_CONTINUE,
        JOB_TYPE,
        JOB_TYPE_NAME,

        JOBOPTION_TYPE,
        JOBOPTION_VARIABLE,
        JOBOPTION_VALUE,
        JOBOPTION_LABEL,
        JOBOPTION_DEFAULT_VALUE,
        JOBOPTION_MINVAL,
        JOBOPTION_MAXVAL,
        JOBOPTION_STEPVAL,
        JOBOPTION_HELPTEXT,
        JOBOPTION_PATTERN,
        JOBOPTION_DIRECTORY,
        JOBOPTION_MENUOPTIONS,

        MATRIX_1_1,
        MATRIX_1_2,
        MATRIX_1_3,
        MATRIX_2_1,
        MATRIX_2_2,
        MATRIX_2_3,
        MATRIX_3_1,
        MATRIX_3_2,
        MATRIX_3_3,

        MICROGRAPH_ACCUM_MOTION_TOTAL,
        MICROGRAPH_ACCUM_MOTION_EARLY,
        MICROGRAPH_ACCUM_MOTION_LATE,
        MICROGRAPH_ID,
        MICROGRAPH_NAME,
        MICROGRAPH_GAIN_NAME,
        MICROGRAPH_DEFECT_FILE,
        MICROGRAPH_NAME_WODOSE,
        MICROGRAPH_MOVIE_NAME,
        MICROGRAPH_METADATA_NAME,
        MICROGRAPH_TILT_ANGLE,
        MICROGRAPH_TILT_AXIS_DIRECTION,
        MICROGRAPH_TILT_AXIS_OUTOFPLANE,
        MICROGRAPH_ORIGINAL_PIXEL_SIZE,
        MICROGRAPH_PIXEL_SIZE,
        MICROGRAPH_PRE_EXPOSURE,
        MICROGRAPH_DOSE_RATE,
        MICROGRAPH_BINNING,
        MICROGRAPH_FRAME_NUMBER,
        MICROGRAPH_MOTION_MODEL_VERSION,
        MICROGRAPH_START_FRAME,
        MICROGRAPH_END_FRAME,
        MICROGRAPH_SHIFT_X,
        MICROGRAPH_SHIFT_Y,
        MICROGRAPH_MOTION_COEFFS_IDX,
        MICROGRAPH_MOTION_COEFF,
        MICROGRAPH_EER_UPSAMPLING,
        MICROGRAPH_EER_GROUPING,

        MASK_NAME,

        MLMODEL_ACCURACY_ROT,
        MLMODEL_ACCURACY_TRANS, // deprecated
        MLMODEL_ACCURACY_TRANS_ANGSTROM,
        MLMODEL_AVE_PMAX,
        MLMODEL_CURRENT_RESOLUTION,
        MLMODEL_CURRENT_SIZE,
        MLMODEL_DATA_VS_PRIOR_REF,
        MLMODEL_DIMENSIONALITY,
        MLMODEL_DIMENSIONALITY_DATA,
        MLMODEL_DIFF2_HALVES_REF,
        MLMODEL_ESTIM_RESOL_REF,
        MLMODEL_FOURIER_COVERAGE_REF,
        MLMODEL_FOURIER_COVERAGE_TOTAL_REF,
        MLMODEL_FSC_HALVES_REF,
        MLMODEL_GROUP_NAME,
        MLMODEL_GROUP_NO,
        MLMODEL_GROUP_NR_PARTICLES,
        MLMODEL_GROUP_SCALE_CORRECTION,
        MLMODEL_HELICAL_NR_ASU,
        MLMODEL_HELICAL_TWIST,
        MLMODEL_HELICAL_TWIST_MIN,
        MLMODEL_HELICAL_TWIST_MAX,
        MLMODEL_HELICAL_TWIST_INITIAL_STEP,
        MLMODEL_HELICAL_RISE,
        MLMODEL_HELICAL_RISE_MIN,
        MLMODEL_HELICAL_RISE_MAX,
        MLMODEL_HELICAL_RISE_INITIAL_STEP,
        MLMODEL_IS_HELIX,
        MLMODEL_INTERPOLATOR,
        MLMODEL_LL,
        MLMODEL_MINIMUM_RADIUS_NN_INTERPOLATION,
        MLMODEL_NORM_CORRECTION_AVG,
        MLMODEL_NR_BODIES,
        MLMODEL_NR_CLASSES,
        MLMODEL_NR_GROUPS,
        MLMODEL_ORIGINAL_SIZE,
        MLMODEL_ORIENTABILITY_CONTRIBUTION,
        MLMODEL_PADDING_FACTOR,
        MLMODEL_PDF_CLASS,
        MLMODEL_PRIOR_OFFX_CLASS,
        MLMODEL_PRIOR_OFFY_CLASS,
        MLMODEL_PDF_ORIENT,
        MLMODEL_PIXEL_SIZE,
        MLMODEL_POWER_REF,
        MLMODEL_PRIOR_MODE,
        MLMODEL_SIGMA_OFFSET, // deprecated
        MLMODEL_SIGMA_OFFSET_ANGSTROM,
        MLMODEL_SIGMA_ROT,
        MLMODEL_SIGMA_TILT,
        MLMODEL_SIGMA_PSI,
        MLMODEL_REF_IMAGE,
        MLMODEL_SGD_GRADIENT_IMAGE,
        MLMODEL_SIGMA2_NOISE,
        MLMODEL_SIGMA2_REF,
        MLMODEL_SSNR_REF,
        MLMODEL_TAU2_FUDGE_FACTOR,
        MLMODEL_TAU2_REF,

        OPTIMISER_ACCURACY_ROT,
        OPTIMISER_ACCURACY_TRANS, // deprecated
        OPTIMISER_ACCURACY_TRANS_ANGSTROM,
        OPTIMISER_ADAPTIVE_FRACTION,
        OPTIMISER_ADAPTIVE_OVERSAMPLING,
        OPTIMISER_AUTO_LOCAL_HP_ORDER,
        OPTIMISER_AVAILABLE_MEMORY,
        OPTIMISER_BEST_RESOL_THUS_FAR,
        OPTIMISER_CHANGES_OPTIMAL_OFFSETS,
        OPTIMISER_CHANGES_OPTIMAL_ORIENTS,
        OPTIMISER_CHANGES_OPTIMAL_CLASSES,
        OPTIMISER_COARSE_SIZE,
        OPTIMISER_DATA_ARE_CTF_PHASE_FLIPPED,
        OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED,
        OPTIMISER_DATA_STARFILE,
        OPTIMISER_DO_AUTO_REFINE,
        OPTIMISER_DO_ONLY_FLIP_CTF_PHASES,
        OPTIMISER_DO_CORRECT_CTF,
        OPTIMISER_DO_CORRECT_MAGNIFICATION,
        OPTIMISER_DO_CORRECT_NORM,
        OPTIMISER_DO_CORRECT_SCALE,
        OPTIMISER_DO_EXTERNAL_RECONSTRUCT,
        OPTIMISER_DO_REALIGN_MOVIES,
        OPTIMISER_DO_MAP,
        OPTIMISER_DO_SGD,
        OPTIMISER_DO_STOCHASTIC_EM,
        OPTIMISER_EXTERNAL_RECONS_DATA_REAL,
        OPTIMISER_EXTERNAL_RECONS_DATA_IMAG,
        OPTIMISER_EXTERNAL_RECONS_WEIGHT,
        OPTIMISER_EXTERNAL_RECONS_RESULT,
        OPTIMISER_EXTERNAL_RECONS_NEWSTAR,
        OPTIMISER_FAST_SUBSETS,
        OPTIMISER_SGD_INI_ITER,
        OPTIMISER_SGD_FIN_ITER,
        OPTIMISER_SGD_INBETWEEN_ITER,
        OPTIMISER_SGD_INI_RESOL,
        OPTIMISER_SGD_FIN_RESOL,
        OPTIMISER_SGD_INI_SUBSET_SIZE,
        OPTIMISER_SGD_FIN_SUBSET_SIZE,
        OPTIMISER_SGD_MU,
        OPTIMISER_SGD_SIGMA2FUDGE_INI,
        OPTIMISER_SGD_SIGMA2FUDGE_HALFLIFE,
        OPTIMISER_SGD_SKIP_ANNNEAL,
        OPTIMISER_SGD_SUBSET_SIZE,
        OPTIMISER_SGD_WRITE_EVERY_SUBSET,
        OPTIMISER_SGD_MAX_SUBSETS,
        OPTIMISER_SGD_STEPSIZE,
        OPTIMISER_DO_SOLVENT_FLATTEN,
        OPTIMISER_DO_SOLVENT_FSC,
        OPTIMISER_DO_SKIP_ALIGN,
        OPTIMISER_DO_SKIP_ROTATE,
        OPTIMISER_DO_SPLIT_RANDOM_HALVES,
        OPTIMISER_DO_ZERO_MASK,
        OPTIMISER_FIX_SIGMA_NOISE,
        OPTIMISER_FIX_SIGMA_OFFSET,
        OPTIMISER_FIX_TAU,
        OPTIMISER_HAS_CONVERGED,
        OPTIMISER_HAS_HIGH_FSC_AT_LIMIT,
        OPTIMISER_HAS_LARGE_INCR_SIZE_ITER_AGO,
        OPTIMISER_DO_HELICAL_REFINE,
        OPTIMISER_IGNORE_HELICAL_SYMMETRY,
        OPTIMISER_FOURIER_MASK,
        OPTIMISER_HELICAL_TWIST_INITIAL,
        OPTIMISER_HELICAL_RISE_INITIAL,
        OPTIMISER_HELICAL_Z_PERCENTAGE,
        OPTIMISER_HELICAL_NSTART,
        OPTIMISER_HELICAL_TUBE_INNER_DIAMETER,
        OPTIMISER_HELICAL_TUBE_OUTER_DIAMETER,
        OPTIMISER_HELICAL_SYMMETRY_LOCAL_REFINEMENT,
        OPTIMISER_HELICAL_SIGMA_DISTANCE,
        OPTIMISER_HELICAL_KEEP_TILT_PRIOR_FIXED,
        OPTIMISER_LOWRES_LIMIT_EXP,
        OPTIMISER_HIGHRES_LIMIT_EXP,
        OPTIMISER_HIGHRES_LIMIT_SGD,
        OPTIMISER_IGNORE_CTF_UNTIL_FIRST_PEAK,
        OPTIMISER_INCR_SIZE,
        OPTIMISER_ITERATION_NO,
        OPTIMISER_LOCAL_SYMMETRY_FILENAME,
        OPTIMISER_LOWRES_JOIN_RANDOM_HALVES,
        OPTIMISER_MAGNIFICATION_RANGE,
        OPTIMISER_MAGNIFICATION_STEP,
        OPTIMISER_MAX_COARSE_SIZE,
        OPTIMISER_MAX_NR_POOL,
        OPTIMISER_MODEL_STARFILE,
        OPTIMISER_MODEL_STARFILE2,
        OPTIMISER_NR_ITERATIONS,
        OPTIMISER_NR_ITER_WO_RESOL_GAIN,
        OPTIMISER_NR_ITER_WO_HIDDEN_VAR_CHANGES,
        OPTIMISER_OPTICS_STARFILE,
        OPTIMISER_OUTPUT_ROOTNAME,
        OPTIMISER_PARTICLE_DIAMETER,
        OPTIMISER_RADIUS_MASK_3D_MAP,
        OPTIMISER_RADIUS_MASK_EXP_PARTICLES,
        OPTIMISER_RANDOM_SEED,
        OPTIMISER_REFS_ARE_CTF_CORRECTED,
        OPTIMISER_SAMPLING_STARFILE,
        OPTIMISER_SMALLEST_CHANGES_OPT_CLASSES,
        OPTIMISER_SMALLEST_CHANGES_OPT_OFFSETS,
        OPTIMISER_SMALLEST_CHANGES_OPT_ORIENTS,
        OPTIMISER_SOLVENT_MASK_NAME,
        OPTIMISER_SOLVENT_MASK2_NAME,
        OPTIMISER_TAU_SPECTRUM_NAME,
        OPTIMISER_USE_TOO_COARSE_SAMPLING,
        OPTIMISER_WIDTH_MASK_EDGE,

        ORIENT_FLIP,
        ORIENT_ID,
        ORIENT_ORIGIN_X, // (deprecated)
        ORIENT_ORIGIN_Y, // (deprecated)
        ORIENT_ORIGIN_Z, // (deprecated)
        ORIENT_ORIGIN_X_PRIOR, // (deprecated)
        ORIENT_ORIGIN_Y_PRIOR, // (deprecated)
        ORIENT_ORIGIN_Z_PRIOR, // (deprecated)
        ORIENT_ORIGIN_X_ANGSTROM,
        ORIENT_ORIGIN_Y_ANGSTROM,
        ORIENT_ORIGIN_Z_ANGSTROM,
        ORIENT_ORIGIN_X_PRIOR_ANGSTROM,
        ORIENT_ORIGIN_Y_PRIOR_ANGSTROM,
        ORIENT_ORIGIN_Z_PRIOR_ANGSTROM,
        ORIENT_ROT,
        ORIENT_ROT_PRIOR,
        ORIENT_ROT_PRIOR_FLIP_RATIO,	// KThurber
        ORIENT_TILT,
        ORIENT_TILT_PRIOR,
        ORIENT_PSI,
        ORIENT_PSI_PRIOR,
        ORIENT_PSI_PRIOR_FLIP_RATIO,
        ORIENT_PSI_PRIOR_FLIP,  // KThurber

        PARTICLE_AUTOPICK_FOM,
        PARTICLE_HELICAL_TUBE_ID,
        PARTICLE_HELICAL_TUBE_PITCH,
        PARTICLE_HELICAL_TRACK_LENGTH, //deprecated
        PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM,
        PARTICLE_CLASS,
        PARTICLE_DLL,
        PARTICLE_ID,
        PARTICLE_FOM,
        PARTICLE_KL_DIVERGENCE,
        PARTICLE_RANDOM_SUBSET,
        PARTICLE_BEAM_TILT_CLASS,
        PARTICLE_NAME,
        PARTICLE_ORI_NAME,
        PARTICLE_NR_SIGNIFICANT_SAMPLES,
        PARTICLE_NR_FRAMES,
        PARTICLE_NR_FRAMES_AVG,
        PARTICLE_MOVIE_RUNNING_AVG,
        PARTICLE_PMAX,
        PARTICLE_NUMBER,

        PIPELINE_JOB_COUNTER,
        PIPELINE_NODE_NAME,
        PIPELINE_NODE_TYPE,
        PIPELINE_PROCESS_ALIAS,
        PIPELINE_PROCESS_NAME,
        PIPELINE_PROCESS_TYPE,
        PIPELINE_PROCESS_STATUS,
        PIPELINE_EDGE_FROM,
        PIPELINE_EDGE_TO,
        PIPELINE_EDGE_PROCESS,

        POSTPROCESS_BFACTOR,
        POSTPROCESS_FINAL_RESOLUTION,
        POSTPROCESS_FRACTION_MOLWEIGHT,
        POSTPROCESS_FRACTION_SOLVENT_MASK,
        POSTPROCESS_FSC_GENERAL,
        POSTPROCESS_FSC_TRUE,
        POSTPROCESS_FSC_PART_MOLWEIGHT,
        POSTPROCESS_FSC_PART_FRACMASK,
        POSTPROCESS_FSC_MASKED,
        POSTPROCESS_FSC_UNMASKED,
        POSTPROCESS_FSC_RANDOM_MASKED,
        POSTPROCESS_AMPLCORR_MASKED,
        POSTPROCESS_AMPLCORR_UNMASKED,
        POSTPROCESS_DPR_MASKED,
        POSTPROCESS_DPR_UNMASKED,
        POSTPROCESS_GUINIER_FIT_CORRELATION,
        POSTPROCESS_GUINIER_FIT_INTERCEPT,
        POSTPROCESS_GUINIER_FIT_SLOPE,
        POSTPROCESS_GUINIER_VALUE_IN,
        POSTPROCESS_GUINIER_VALUE_INVMTF,
        POSTPROCESS_GUINIER_VALUE_WEIGHTED,
        POSTPROCESS_GUINIER_VALUE_SHARPENED,
        POSTPROCESS_GUINIER_VALUE_INTERCEPT,
        POSTPROCESS_GUINIER_RESOL_SQUARED,
        POSTPROCESS_MOLWEIGHT,
        POSTPROCESS_MTF_VALUE, ///< Detector MTF value
        POSTPROCESS_RANDOMISE_FROM,
        POSTPROCESS_UNFIL_HALFMAP1,
        POSTPROCESS_UNFIL_HALFMAP2,

        SAMPLING_IS_3D,
        SAMPLING_IS_3D_TRANS,
        SAMPLING_HEALPIX_ORDER,
        SAMPLING_HEALPIX_ORDER_ORI,
        SAMPLING_LIMIT_TILT,
        SAMPLING_OFFSET_RANGE,
        SAMPLING_OFFSET_STEP,
        SAMPLING_OFFSET_RANGE_ORI,
        SAMPLING_OFFSET_STEP_ORI,
        SAMPLING_HELICAL_OFFSET_STEP,
        SAMPLING_PERTURB,
        SAMPLING_PERTURBATION_FACTOR,
        SAMPLING_PRIOR_MODE,
        SAMPLING_PSI_STEP,
        SAMPLING_PSI_STEP_ORI,
        SAMPLING_SIGMA_ROT,
        SAMPLING_SIGMA_TILT,
        SAMPLING_SIGMA_PSI,
        SAMPLING_SYMMETRY,

        SCHEDULE_EDGE_NUMBER,
        SCHEDULE_EDGE_INPUT,
        SCHEDULE_EDGE_OUTPUT,
        SCHEDULE_EDGE_IS_FORK,
        SCHEDULE_EDGE_OUTPUT_TRUE,
        SCHEDULE_EDGE_BOOLEAN,
        SCHEDULE_GENERAL_CURRENT_NODE,
        SCHEDULE_GENERAL_ORIGINAL_START_NODE,
        SCHEDULE_GENERAL_EMAIL,
        SCHEDULE_GENERAL_NAME,
        SCHEDULE_JOB_NAME,
        SCHEDULE_JOB_ORI_NAME,
        SCHEDULE_JOB_MODE,
        SCHEDULE_JOB_HAS_STARTED,
        SCHEDULE_OPERATOR_NAME,
        SCHEDULE_OPERATOR_TYPE,
        SCHEDULE_OPERATOR_INPUT1,
        SCHEDULE_OPERATOR_INPUT2,
        SCHEDULE_OPERATOR_OUTPUT,
        SCHEDULE_VAR_BOOL_NAME,
        SCHEDULE_VAR_BOOL_VALUE,
        SCHEDULE_VAR_BOOL_ORI_VALUE,
        SCHEDULE_VAR_FLOAT_NAME,
        SCHEDULE_VAR_FLOAT_VALUE,
        SCHEDULE_VAR_FLOAT_ORI_VALUE,
        SCHEDULE_VAR_STRING_NAME,
        SCHEDULE_VAR_STRING_VALUE,
        SCHEDULE_VAR_STRING_ORI_VALUE,

        SELECTED,
        SELECT_PARTICLES_ZSCORE,
        SORTED_IDX,
        STARFILE_MOVIE_PARTICLES,
        PERFRAME_CUMULATIVE_WEIGHT,
        PERFRAME_RELATIVE_WEIGHT,

        RESOLUTION,
        RESOLUTION_ANGSTROM,
        RESOLUTION_INVPIXEL,
        SPECTRAL_IDX,

        UNKNOWN_LABEL,

        LAST_LABEL
        /** NOTE: Keep this label at the end.
         * It is here for looping purposes.
         */
    };

    static EMDLabel str2Label(const std::string &labelName);
    static std::string label2Str(const EMDLabel &label);

    template <typename T>
    static bool is(const EMDLabel &label);

    static bool isValidLabel(const EMDLabel &label);
    static bool isValidLabel(const std::string &labelName);

    static void printDefinitions(std::ostream &out);

    private:

    static std::map<const EMDLabel, const LabelData> data;
    static std::map<std::string, EMDLabel> labels;
    static std::map<std::string, std::string> definitions;
    static StaticInitialization initialization; // Just for initialization

    template <typename T>
    static void addLabel(EMDLabel label, const std::string &name, const std::string &definition = "undocumented");
    static void addAltLabel(EMDLabel label, std::string name);

    friend class StaticInitialization;

};

// Just a class for static initialization
class StaticInitialization {

    private:

    /** Call addLabel on every value in the EMDLabel enum,
     *  giving each a type, a name, and a definition.
     */
    StaticInitialization() {
        /// NOTE: ==== Add labels entries from here in the SAME ORDER as declared in ENUM ==========
        EMDL::addLabel<std::string>(EMDL::COMMENT, "rlnComment", "A metadata comment (This is treated in a special way)");

        EMDL::addLabel<int>        (EMDL::AREA_ID,   "rlnAreaId", "ID (i.e. a unique number) of an area (i.e. field-of-view)");
        EMDL::addLabel<std::string>(EMDL::AREA_NAME, "rlnAreaName", "Name of an area (i.e. field-of-view)");

        EMDL::addLabel<std::string>(EMDL::BODY_MASK_NAME, "rlnBodyMaskName", "Name of an image that contains a [0,1] body mask for multi-body refinement");
        EMDL::addLabel<int>(EMDL::BODY_KEEP_FIXED,    "rlnBodyKeepFixed", "Flag to indicate whether to keep a body fixed (value 1) or keep on refining it (0)");
        EMDL::addLabel<std::string>(EMDL::BODY_REFERENCE_NAME, "rlnBodyReferenceName", "Name of an image that contains the initial reference for one body of a multi-body refinement");
        EMDL::addLabel<double>(EMDL::BODY_ROTATE_DIRECTION_X, "rlnBodyRotateDirectionX", "X-component of axis around which to rotate this body");
        EMDL::addLabel<double>(EMDL::BODY_ROTATE_DIRECTION_Y, "rlnBodyRotateDirectionY", "Y-component of axis around which to rotate this body");
        EMDL::addLabel<double>(EMDL::BODY_ROTATE_DIRECTION_Z, "rlnBodyRotateDirectionZ", "Z-component of axis around which to rotate this body");
        EMDL::addLabel<int>(EMDL::BODY_ROTATE_RELATIVE_TO,    "rlnBodyRotateRelativeTo", "Number of the body relative to which this body rotates (if negative, use rlnBodyRotateDirectionXYZ)");
        EMDL::addLabel<double>(EMDL::BODY_SIGMA_ANG, "rlnBodySigmaAngles", "Width of prior on all three Euler angles of a body in multibody refinement (in degrees)");
        EMDL::addLabel<double>(EMDL::BODY_SIGMA_OFFSET, "rlnBodySigmaOffset", "Width of prior on origin offsets of a body in multibody refinement (in pixels)");
        EMDL::addLabel<double>(EMDL::BODY_SIGMA_OFFSET_ANGSTROM, "rlnBodySigmaOffsetAngst", "Width of prior on origin offsets of a body in multibody refinement (in Angstroms)");
        EMDL::addLabel<double>(EMDL::BODY_SIGMA_ROT, "rlnBodySigmaRot", "Width of prior on rot angles of a body in multibody refinement (in degrees)");
        EMDL::addLabel<double>(EMDL::BODY_SIGMA_TILT, "rlnBodySigmaTilt", "Width of prior on tilt angles of a body in multibody refinement (in degrees)");
        EMDL::addLabel<double>(EMDL::BODY_SIGMA_PSI, "rlnBodySigmaPsi", "Width of prior on psi angles of a body in multibody refinement (in degrees)");
        EMDL::addLabel<std::string>(EMDL::BODY_STAR_FILE, "rlnBodyStarFile", "Name of STAR file with body masks and metadata");

        EMDL::addLabel<double>(EMDL::CTF_ASTIGMATISM, "rlnCtfAstigmatism", "Absolute value of the difference between defocus in U- and V-direction (in A)");
        EMDL::addLabel<double>(EMDL::CTF_BFACTOR, "rlnCtfBfactor", "B-factor (in A^2) that describes CTF power spectrum fall-off");
        EMDL::addLabel<double>(EMDL::CTF_MAXRES, "rlnCtfMaxResolution", "Estimated maximum resolution (in A) of significant CTF Thon rings");
        EMDL::addLabel<double>(EMDL::CTF_VALIDATIONSCORE, "rlnCtfValidationScore", "Gctf-based validation score for the quality of the CTF fit");
        EMDL::addLabel<double>(EMDL::CTF_SCALEFACTOR, "rlnCtfScalefactor", "Linear scale-factor on the CTF (values between 0 and 1)");
        EMDL::addLabel<double>(EMDL::CTF_VOLTAGE, "rlnVoltage", "Voltage of the microscope (in kV)");
        EMDL::addLabel<double>(EMDL::CTF_DEFOCUSU, "rlnDefocusU", "Defocus in U-direction (in Angstroms, positive values for underfocus)");
        EMDL::addLabel<double>(EMDL::CTF_DEFOCUSV, "rlnDefocusV", "Defocus in V-direction (in Angstroms, positive values for underfocus)");
        EMDL::addLabel<double>(EMDL::CTF_DEFOCUS_ANGLE, "rlnDefocusAngle", "Angle between X and defocus U direction (in degrees)");
        EMDL::addLabel<double>(EMDL::CTF_CS, "rlnSphericalAberration", "Spherical aberration (in millimeters)");
        EMDL::addLabel<double>(EMDL::CTF_CA, "rlnChromaticAberration", "Chromatic aberration (in millimeters)");
        EMDL::addLabel<double>(EMDL::CTF_DETECTOR_PIXEL_SIZE, "rlnDetectorPixelSize", "Pixel size of the detector (in micrometers)");
        EMDL::addLabel<std::string>(EMDL::CTF_POWER_SPECTRUM, "rlnCtfPowerSpectrum", "Power spectrum for CTF estimation");
        EMDL::addLabel<double>(EMDL::CTF_ENERGY_LOSS, "rlnEnergyLoss", "Energy loss (in eV)");
        EMDL::addLabel<double>(EMDL::CTF_FOM, "rlnCtfFigureOfMerit", "Figure of merit for the fit of the CTF (not used inside relion_refine)");
        EMDL::addLabel<std::string>(EMDL::CTF_IMAGE, "rlnCtfImage", "Name of an image with all CTF values");
        EMDL::addLabel<double>(EMDL::CTF_LENS_STABILITY, "rlnLensStability", "Lens stability (in ppm)");
        EMDL::addLabel<double>(EMDL::CTF_MAGNIFICATION, "rlnMagnification", "Magnification at the detector (in times)");
        EMDL::addLabel<double>(EMDL::CTF_PHASESHIFT, "rlnPhaseShift", "Phase-shift from a phase-plate (in degrees)");
        EMDL::addLabel<double>(EMDL::CTF_CONVERGENCE_CONE, "rlnConvergenceCone", "Convergence cone (in mrad)");
        EMDL::addLabel<double>(EMDL::CTF_LONGITUDINAL_DISPLACEMENT, "rlnLongitudinalDisplacement", "Longitudinal displacement (in Angstroms)");
        EMDL::addLabel<double>(EMDL::CTF_TRANSVERSAL_DISPLACEMENT, "rlnTransversalDisplacement", "Transversal displacement (in Angstroms)");
        EMDL::addLabel<double>(EMDL::CTF_Q0, "rlnAmplitudeContrast", "Amplitude contrast (as a fraction, i.e. 10% = 0.1)");
        EMDL::addLabel<double>(EMDL::CTF_VALUE, "rlnCtfValue", "Value of the Contrast Transfer Function");

        EMDL::addLabel<std::string>(EMDL::IMAGE_NAME,        "rlnImageName", "Name of an image");
        EMDL::addLabel<std::string>(EMDL::IMAGE_ORI_NAME,        "rlnImageOriginalName", "Original name of an image");
        EMDL::addLabel<std::string>(EMDL::IMAGE_RECONSTRUCT_NAME,        "rlnReconstructImageName", "Name of an image to be used for reconstruction only");
        EMDL::addLabel<int>(EMDL::IMAGE_ID,           "rlnImageId", "ID (i.e. a unique number) of an image");
        EMDL::addLabel<bool>(EMDL::IMAGE_ENABLED,          "rlnEnabled", "Not used in RELION, only included for backward compatibility with XMIPP selfiles");
        EMDL::addLabel<int>(EMDL::IMAGE_DATATYPE,           "rlnDataType", "Type of data stored in an image (e.g. int, RFLOAT etc)");
        EMDL::addLabel<int>(EMDL::IMAGE_DIMENSIONALITY,           "rlnImageDimensionality", "Dimensionality of data stored in an image (i.e. 2 or 3)");
        EMDL::addLabel<double>(EMDL::IMAGE_BEAMTILT_X,        "rlnBeamTiltX", "Beam tilt in the X-direction (in mrad)");
        EMDL::addLabel<double>(EMDL::IMAGE_BEAMTILT_Y,        "rlnBeamTiltY", "Beam tilt in the Y-direction (in mrad)");
        EMDL::addLabel<std::string>(EMDL::IMAGE_MTF_FILENAME,        "rlnMtfFileName", "The filename of a STAR file with the MTF for this optics group or image");
        EMDL::addLabel<int>(EMDL::IMAGE_OPTICS_GROUP,           "rlnOpticsGroup", "Group of particles with identical optical properties");
        EMDL::addLabel<std::string>(EMDL::IMAGE_OPTICS_GROUP_NAME,        "rlnOpticsGroupName", "The name of a group of particles with identical optical properties");
        EMDL::addLabel<std::vector<double> >(EMDL::IMAGE_ODD_ZERNIKE_COEFFS, "rlnOddZernike", "Coefficients for the antisymmetrical Zernike polynomials");
        EMDL::addLabel<std::vector<double> >(EMDL::IMAGE_EVEN_ZERNIKE_COEFFS, "rlnEvenZernike", "Coefficients for the symmetrical Zernike polynomials");
        EMDL::addLabel<double>(EMDL::IMAGE_PIXEL_SIZE,        "rlnImagePixelSize", "Pixel size (in Angstrom)");
        EMDL::addLabel<double>(EMDL::IMAGE_MAG_MATRIX_00,        "rlnMagMat00", "Anisotropic magnification matrix, element 1,1");
        EMDL::addLabel<double>(EMDL::IMAGE_MAG_MATRIX_01,        "rlnMagMat01", "Anisotropic magnification matrix, element 1,2");
        EMDL::addLabel<double>(EMDL::IMAGE_MAG_MATRIX_10,        "rlnMagMat10", "Anisotropic magnification matrix, element 2,1");
        EMDL::addLabel<double>(EMDL::IMAGE_MAG_MATRIX_11,        "rlnMagMat11", "Anisotropic magnification matrix, element 2,2");

        EMDL::addLabel<double>(EMDL::IMAGE_COORD_X, "rlnCoordinateX", "X-Position of an image in a micrograph (in pixels)");
        EMDL::addLabel<double>(EMDL::IMAGE_COORD_Y, "rlnCoordinateY", "Y-Position of an image in a micrograph (in pixels)");
        EMDL::addLabel<double>(EMDL::IMAGE_COORD_Z, "rlnCoordinateZ", "Z-Position of an image in a 3D micrograph, i.e. tomogram (in pixels)");
        EMDL::addLabel<int>(EMDL::IMAGE_FRAME_NR,    "rlnMovieFrameNumber", "Number of a movie frame");
        EMDL::addLabel<double>(EMDL::IMAGE_NORM_CORRECTION, "rlnNormCorrection", "Normalisation correction value for an image");
        EMDL::addLabel<double>(EMDL::IMAGE_MAGNIFICATION_CORRECTION, "rlnMagnificationCorrection", "Magnification correction value for an image");
        EMDL::addLabel<double>(EMDL::IMAGE_SAMPLINGRATE, "rlnSamplingRate", "Sampling rate of an image (in Angstrom/pixel)");
        EMDL::addLabel<double>(EMDL::IMAGE_SAMPLINGRATE_X, "rlnSamplingRateX", "Sampling rate in X-direction of an image (in Angstrom/pixel)");
        EMDL::addLabel<double>(EMDL::IMAGE_SAMPLINGRATE_Y, "rlnSamplingRateY", "Sampling rate in Y-direction of an image (in Angstrom/pixel)");
        EMDL::addLabel<double>(EMDL::IMAGE_SAMPLINGRATE_Z, "rlnSamplingRateZ", "Sampling rate in Z-direction of an image (in Angstrom/pixel)");
        EMDL::addLabel<int>(EMDL::IMAGE_SIZE,    "rlnImageSize", "Size of an image (in pixels)");
        EMDL::addLabel<int>(EMDL::IMAGE_SIZE_X,    "rlnImageSizeX", "Size of an image in the X-direction (in pixels)");
        EMDL::addLabel<int>(EMDL::IMAGE_SIZE_Y,    "rlnImageSizeY", "Size of an image in the Y-direction (in pixels)");
        EMDL::addLabel<int>(EMDL::IMAGE_SIZE_Z,    "rlnImageSizeZ", "Size of an image in the Z-direction (in pixels)");
        EMDL::addLabel<double>(EMDL::IMAGE_STATS_MIN, "rlnMinimumValue", "Minimum value for the pixels in an image");
        EMDL::addLabel<double>(EMDL::IMAGE_STATS_MAX, "rlnMaximumValue", "Maximum value for the pixels in an image");
        EMDL::addLabel<double>(EMDL::IMAGE_STATS_AVG, "rlnAverageValue", "Average value for the pixels in an image");
        EMDL::addLabel<double>(EMDL::IMAGE_STATS_STDDEV, "rlnStandardDeviationValue", "Standard deviation for the pixel values in an image");
        EMDL::addLabel<double>(EMDL::IMAGE_STATS_SKEW, "rlnSkewnessValue", "Skewness (3rd moment) for the pixel values in an image");
        EMDL::addLabel<double>(EMDL::IMAGE_STATS_KURT, "rlnKurtosisExcessValue", "Kurtosis excess (4th moment - 3) for the pixel values in an image");
        EMDL::addLabel<double>(EMDL::IMAGE_WEIGHT, "rlnImageWeight", "Relative weight of an image");

        EMDL::addLabel<std::string>(EMDL::MASK_NAME, "rlnMaskName", "Name of an image that contains a [0,1] mask");

        EMDL::addLabel<bool>(EMDL::JOB_IS_CONTINUE,   "rlnJobIsContinue", "Is tthis a continuation job?");
        EMDL::addLabel<int>(EMDL::JOB_TYPE,    "rlnJobType", "Which type of job is this?");
        EMDL::addLabel<std::string>(EMDL::JOB_TYPE_NAME, "rlnJobTypeName", "The name for this type of job (also name of main directory for output jobs)");

        EMDL::addLabel<int>(EMDL::JOBOPTION_TYPE, "rlnJoboptionType", "Which type of joboption is this?");
        EMDL::addLabel<std::string>(EMDL::JOBOPTION_VARIABLE, "rlnJobOptionVariable", "Name of the joboption variable");
        EMDL::addLabel<std::string>(EMDL::JOBOPTION_VALUE, "rlnJobOptionValue", "Value of a joboption");
        EMDL::addLabel<std::string>(EMDL::JOBOPTION_LABEL, "rlnJobOptionGUILabel", "GUI label of a joboption");
        EMDL::addLabel<std::string>(EMDL::JOBOPTION_DEFAULT_VALUE, "rlnJobOptionDefaultValue", "Default value of a joboption");
        EMDL::addLabel<double>(EMDL::JOBOPTION_MINVAL, "rlnJobOptionSliderMin", "Minimum value for slider of a joboption");
        EMDL::addLabel<double>(EMDL::JOBOPTION_MAXVAL, "rlnJobOptionSliderMax", "Maximum value for slider of a joboption");
        EMDL::addLabel<double>(EMDL::JOBOPTION_STEPVAL, "rlnJobOptionSliderStep", "Step value for slider of a joboption");
        EMDL::addLabel<std::string>(EMDL::JOBOPTION_HELPTEXT, "rlnJobOptionHelpText", "Extra helptext of a joboption");
        EMDL::addLabel<std::string>(EMDL::JOBOPTION_PATTERN, "rlnJobOptionFilePattern", "Pattern for file browser of a joboption");
        EMDL::addLabel<std::string>(EMDL::JOBOPTION_DIRECTORY, "rlnJobOptionDirectoryDefault", "Default directory for file browser of a joboption");
        EMDL::addLabel<std::string>(EMDL::JOBOPTION_MENUOPTIONS, "rlnJobOptionMenuOptions", "Options for pull-down menu");

        EMDL::addLabel<double>(EMDL::MATRIX_1_1, "rlnMatrix_1_1", "Matrix element (1,1) of a 3x3 matrix");
        EMDL::addLabel<double>(EMDL::MATRIX_1_2, "rlnMatrix_1_2", "Matrix element (1,2) of a 3x3 matrix");
        EMDL::addLabel<double>(EMDL::MATRIX_1_3, "rlnMatrix_1_3", "Matrix element (1,3) of a 3x3 matrix");
        EMDL::addLabel<double>(EMDL::MATRIX_2_1, "rlnMatrix_2_1", "Matrix element (2,1) of a 3x3 matrix");
        EMDL::addLabel<double>(EMDL::MATRIX_2_2, "rlnMatrix_2_2", "Matrix element (2,1) of a 3x3 matrix");
        EMDL::addLabel<double>(EMDL::MATRIX_2_3, "rlnMatrix_2_3", "Matrix element (2,1) of a 3x3 matrix");
        EMDL::addLabel<double>(EMDL::MATRIX_3_1, "rlnMatrix_3_1", "Matrix element (3,1) of a 3x3 matrix");
        EMDL::addLabel<double>(EMDL::MATRIX_3_2, "rlnMatrix_3_2", "Matrix element (3,1) of a 3x3 matrix");
        EMDL::addLabel<double>(EMDL::MATRIX_3_3, "rlnMatrix_3_3", "Matrix element (3,1) of a 3x3 matrix");

        EMDL::addLabel<double>(EMDL::MICROGRAPH_ACCUM_MOTION_TOTAL, "rlnAccumMotionTotal", "Accumulated global motion during the entire movie (in A)");
        EMDL::addLabel<double>(EMDL::MICROGRAPH_ACCUM_MOTION_EARLY, "rlnAccumMotionEarly", "Accumulated global motion during the first frames of the movie (in A)");
        EMDL::addLabel<double>(EMDL::MICROGRAPH_ACCUM_MOTION_LATE, "rlnAccumMotionLate", "Accumulated global motion during the last frames of the movie (in A)");
        EMDL::addLabel<int>(EMDL::MICROGRAPH_ID,    "rlnMicrographId", "ID (i.e. a unique number) of a micrograph");
        EMDL::addLabel<std::string>(EMDL::MICROGRAPH_NAME, "rlnMicrographName", "Name of a micrograph");
        EMDL::addLabel<std::string>(EMDL::MICROGRAPH_GAIN_NAME, "rlnMicrographGainName", "Name of a gain reference");
        EMDL::addLabel<std::string>(EMDL::MICROGRAPH_DEFECT_FILE, "rlnMicrographDefectFile", "Name of a defect list file");
        EMDL::addLabel<std::string>(EMDL::MICROGRAPH_NAME_WODOSE, "rlnMicrographNameNoDW", "Name of a micrograph without dose weighting");
        EMDL::addLabel<std::string>(EMDL::MICROGRAPH_MOVIE_NAME, "rlnMicrographMovieName", "Name of a micrograph movie stack");
        EMDL::addLabel<std::string>(EMDL::MICROGRAPH_METADATA_NAME, "rlnMicrographMetadata", "Name of a micrograph metadata file");
        EMDL::addLabel<double>(EMDL::MICROGRAPH_TILT_ANGLE, "rlnMicrographTiltAngle", "Tilt angle (in degrees) used to collect a micrograph");
        EMDL::addLabel<double>(EMDL::MICROGRAPH_TILT_AXIS_DIRECTION, "rlnMicrographTiltAxisDirection", "Direction of the tilt-axis (in degrees) used to collect a micrograph");
        EMDL::addLabel<double>(EMDL::MICROGRAPH_TILT_AXIS_OUTOFPLANE, "rlnMicrographTiltAxisOutOfPlane", "Out-of-plane angle (in degrees) of the tilt-axis used to collect a micrograph (90=in-plane)");
        EMDL::addLabel<double>(EMDL::MICROGRAPH_ORIGINAL_PIXEL_SIZE, "rlnMicrographOriginalPixelSize", "Pixel size of original movie before binning in Angstrom/pixel.");
        EMDL::addLabel<double>(EMDL::MICROGRAPH_PIXEL_SIZE, "rlnMicrographPixelSize", "Pixel size of (averaged) micrographs after binning in Angstrom/pixel.");
        EMDL::addLabel<double>(EMDL::MICROGRAPH_PRE_EXPOSURE, "rlnMicrographPreExposure", "Pre-exposure dose in electrons per square Angstrom");
        EMDL::addLabel<double>(EMDL::MICROGRAPH_DOSE_RATE, "rlnMicrographDoseRate", "Dose rate in electrons per square Angstrom per frame");
        EMDL::addLabel<double>(EMDL::MICROGRAPH_BINNING, "rlnMicrographBinning", "Micrograph binning factor");
        EMDL::addLabel<int>(EMDL::MICROGRAPH_FRAME_NUMBER,    "rlnMicrographFrameNumber", "Micrograph frame number");
        EMDL::addLabel<int>(EMDL::MICROGRAPH_MOTION_MODEL_VERSION,    "rlnMotionModelVersion", "Version of micrograph motion model");
        EMDL::addLabel<int>(EMDL::MICROGRAPH_START_FRAME,    "rlnMicrographStartFrame", "Start frame of a motion model");
        EMDL::addLabel<int>(EMDL::MICROGRAPH_END_FRAME,    "rlnMicrographEndFrame", "End frame of a motion model");
        EMDL::addLabel<double>(EMDL::MICROGRAPH_SHIFT_X, "rlnMicrographShiftX", "X shift of a (patch of) micrograph");
        EMDL::addLabel<double>(EMDL::MICROGRAPH_SHIFT_Y, "rlnMicrographShiftY", "Y shift of a (patch of) micrograph");
        EMDL::addLabel<int>(EMDL::MICROGRAPH_MOTION_COEFFS_IDX,    "rlnMotionModelCoeffsIdx", "Index of a coefficient of a motion model");
        EMDL::addLabel<double>(EMDL::MICROGRAPH_MOTION_COEFF, "rlnMotionModelCoeff", "A coefficient of a motion model");
        EMDL::addLabel<int>(EMDL::MICROGRAPH_EER_UPSAMPLING,    "rlnEERUpsampling", "EER upsampling ratio (1 = 4K, 2 = 8K)");
        EMDL::addLabel<int>(EMDL::MICROGRAPH_EER_GROUPING,    "rlnEERGrouping", "The number of hardware frames to group");

        EMDL::addLabel<double>(EMDL::MLMODEL_ACCURACY_ROT, "rlnAccuracyRotations", "Estimated accuracy (in degrees) with which rotations can be assigned");
        EMDL::addLabel<double>(EMDL::MLMODEL_ACCURACY_TRANS, "rlnAccuracyTranslations", "Estimated accuracy (in pixels) with which translations can be assigned");
        EMDL::addLabel<double>(EMDL::MLMODEL_ACCURACY_TRANS_ANGSTROM, "rlnAccuracyTranslationsAngst", "Estimated accuracy (in Angstroms) with which translations can be assigned");
        EMDL::addLabel<double>(EMDL::MLMODEL_AVE_PMAX, "rlnAveragePmax", "Average value (over all images) of the maxima of the probability distributions");
        EMDL::addLabel<double>(EMDL::MLMODEL_CURRENT_RESOLUTION, "rlnCurrentResolution", "Current resolution where SSNR^MAP drops below 1 (in 1/Angstroms)");
        EMDL::addLabel<int>(EMDL::MLMODEL_CURRENT_SIZE, "rlnCurrentImageSize", "Current size of the images used in the refinement");
        EMDL::addLabel<double>(EMDL::MLMODEL_DATA_VS_PRIOR_REF, "rlnSsnrMap", "Spectral signal-to-noise ratio as defined for MAP estimation (SSNR^MAP)");
        EMDL::addLabel<int>(EMDL::MLMODEL_DIMENSIONALITY, "rlnReferenceDimensionality", "Dimensionality of the references (2D/3D)");
        EMDL::addLabel<int>(EMDL::MLMODEL_DIMENSIONALITY_DATA, "rlnDataDimensionality", "Dimensionality of the data (2D/3D)");
        EMDL::addLabel<double>(EMDL::MLMODEL_DIFF2_HALVES_REF, "rlnDiff2RandomHalves", "Power of the differences between two independent reconstructions from random halves of the data");
        EMDL::addLabel<double>(EMDL::MLMODEL_ESTIM_RESOL_REF, "rlnEstimatedResolution", "Estimated resolution (in A) for a reference");
        EMDL::addLabel<double>(EMDL::MLMODEL_FOURIER_COVERAGE_REF, "rlnFourierCompleteness", "Fraction of Fourier components (per resolution shell) with SNR>1");
        EMDL::addLabel<double>(EMDL::MLMODEL_FOURIER_COVERAGE_TOTAL_REF, "rlnOverallFourierCompleteness", "Fraction of all Fourier components up to the current resolution with SNR>1");
        EMDL::addLabel<double>(EMDL::MLMODEL_FSC_HALVES_REF, "rlnGoldStandardFsc", "Fourier shell correlation between two independent reconstructions from random halves of the data");
        EMDL::addLabel<std::string>(EMDL::MLMODEL_GROUP_NAME, "rlnGroupName", "The name of a group of images (e.g. all images from a micrograph)");
        EMDL::addLabel<int>(EMDL::MLMODEL_GROUP_NO, "rlnGroupNumber", "The number of a group of images");
        EMDL::addLabel<int>(EMDL::MLMODEL_GROUP_NR_PARTICLES, "rlnGroupNrParticles", "Number particles in a group of images");
        EMDL::addLabel<double>(EMDL::MLMODEL_GROUP_SCALE_CORRECTION, "rlnGroupScaleCorrection", "Intensity-scale correction for a group of images");
        EMDL::addLabel<int>(EMDL::MLMODEL_HELICAL_NR_ASU, "rlnNrHelicalAsymUnits", "How many new helical asymmetric units are there in each box");
        EMDL::addLabel<double>(EMDL::MLMODEL_HELICAL_TWIST, "rlnHelicalTwist", "The helical twist (rotation per subunit) in degrees");
        EMDL::addLabel<double>(EMDL::MLMODEL_HELICAL_TWIST_MIN, "rlnHelicalTwistMin", "Minimum helical twist (in degrees, + for right-handedness)");
        EMDL::addLabel<double>(EMDL::MLMODEL_HELICAL_TWIST_MAX, "rlnHelicalTwistMax", "Maximum helical twist (in degrees, + for right-handedness)");
        EMDL::addLabel<double>(EMDL::MLMODEL_HELICAL_TWIST_INITIAL_STEP, "rlnHelicalTwistInitialStep", "Initial step of helical twist search (in degrees)");
        EMDL::addLabel<double>(EMDL::MLMODEL_HELICAL_RISE, "rlnHelicalRise", "The helical rise (translation per subunit) in Angstroms");
        EMDL::addLabel<double>(EMDL::MLMODEL_HELICAL_RISE_MIN, "rlnHelicalRiseMin", "Minimum helical rise (in Angstroms)");
        EMDL::addLabel<double>(EMDL::MLMODEL_HELICAL_RISE_MAX, "rlnHelicalRiseMax", "Maximum helical rise (in Angstroms)");
        EMDL::addLabel<double>(EMDL::MLMODEL_HELICAL_RISE_INITIAL_STEP, "rlnHelicalRiseInitialStep", "Initial step of helical rise search (in Angstroms)");
        EMDL::addLabel<bool>(EMDL::MLMODEL_IS_HELIX, "rlnIsHelix", "Flag to indicate that helical refinement should be performed");
        EMDL::addLabel<int>(EMDL::MLMODEL_INTERPOLATOR, "rlnFourierSpaceInterpolator", "The kernel used for Fourier-space interpolation (NN=0, linear=1)");
        EMDL::addLabel<double>(EMDL::MLMODEL_LL, "rlnLogLikelihood", "Value of the log-likelihood target function");
        EMDL::addLabel<int>(EMDL::MLMODEL_MINIMUM_RADIUS_NN_INTERPOLATION, "rlnMinRadiusNnInterpolation", "Minimum radius for NN-interpolation (in Fourier pixels), for smaller radii linear int. is used");
        EMDL::addLabel<double>(EMDL::MLMODEL_NORM_CORRECTION_AVG, "rlnNormCorrectionAverage", "Average value (over all images) of the normalisation correction values");
        EMDL::addLabel<int>(EMDL::MLMODEL_NR_CLASSES, "rlnNrClasses", "The number of references (i.e. classes) to be used in refinement");
        EMDL::addLabel<int>(EMDL::MLMODEL_NR_BODIES, "rlnNrBodies", "The number of independent rigid bodies to be refined in multi-body refinement");
        EMDL::addLabel<int>(EMDL::MLMODEL_NR_GROUPS, "rlnNrGroups", "The number of different groups of images (each group has its own noise spectrum, and intensity-scale correction)");
        EMDL::addLabel<double>(EMDL::MLMODEL_ORIENTABILITY_CONTRIBUTION, "rlnSpectralOrientabilityContribution", "Spectral SNR contribution to the orientability of individual particles");
        EMDL::addLabel<int>(EMDL::MLMODEL_ORIGINAL_SIZE, "rlnOriginalImageSize", "Original size of the images (in pixels)");
        EMDL::addLabel<double>(EMDL::MLMODEL_PADDING_FACTOR, "rlnPaddingFactor", "Oversampling factor for Fourier transforms of the references");
        EMDL::addLabel<double>(EMDL::MLMODEL_PDF_CLASS, "rlnClassDistribution", "Probability Density Function of the different classes (i.e. fraction of images assigned to each class)");
        EMDL::addLabel<double>(EMDL::MLMODEL_PRIOR_OFFX_CLASS, "rlnClassPriorOffsetX", "Prior in the X-offset for a class (in pixels)");
        EMDL::addLabel<double>(EMDL::MLMODEL_PRIOR_OFFY_CLASS, "rlnClassPriorOffsetY", "Prior in the Y-offset for a class (in pixels)");
        EMDL::addLabel<double>(EMDL::MLMODEL_PDF_ORIENT, "rlnOrientationDistribution", "Probability Density Function of the orientations  (i.e. fraction of images assigned to each orient)");
        EMDL::addLabel<double>(EMDL::MLMODEL_PIXEL_SIZE, "rlnPixelSize", "Size of the pixels in the references and images (in Angstroms)");
        EMDL::addLabel<double>(EMDL::MLMODEL_POWER_REF, "rlnReferenceSpectralPower", "Spherical average of the power of the reference");
        EMDL::addLabel<int>(EMDL::MLMODEL_PRIOR_MODE, "rlnOrientationalPriorMode", "Mode for prior distributions on the orientations (0=no prior; 1=(rot,tilt,psi); 2=(rot,tilt); 3=rot; 4=tilt; 5=psi) ");
        EMDL::addLabel<std::string>(EMDL::MLMODEL_REF_IMAGE, "rlnReferenceImage", "Name of a reference image");
        EMDL::addLabel<std::string>(EMDL::MLMODEL_SGD_GRADIENT_IMAGE, "rlnSGDGradientImage", "Name of image containing the SGD gradient");
        EMDL::addLabel<double>(EMDL::MLMODEL_SIGMA_OFFSET, "rlnSigmaOffsets", "Standard deviation in the origin offsets (in pixels)");
        EMDL::addLabel<double>(EMDL::MLMODEL_SIGMA_OFFSET_ANGSTROM, "rlnSigmaOffsetsAngst", "Standard deviation in the origin offsets (in Angstroms)");
        EMDL::addLabel<double>(EMDL::MLMODEL_SIGMA2_NOISE, "rlnSigma2Noise", "Spherical average of the standard deviation in the noise (sigma)");
        EMDL::addLabel<double>(EMDL::MLMODEL_SIGMA2_REF, "rlnReferenceSigma2", "Spherical average of the estimated power in the noise of a reference");
        EMDL::addLabel<double>(EMDL::MLMODEL_SIGMA_ROT, "rlnSigmaPriorRotAngle", "Standard deviation of the prior on the rot (i.e. first Euler) angle");
        EMDL::addLabel<double>(EMDL::MLMODEL_SIGMA_TILT, "rlnSigmaPriorTiltAngle", "Standard deviation of the prior on the tilt (i.e. second Euler) angle");
        EMDL::addLabel<double>(EMDL::MLMODEL_SIGMA_PSI, "rlnSigmaPriorPsiAngle", "Standard deviation of the prior on the psi (i.e. third Euler) angle");
        EMDL::addLabel<double>(EMDL::MLMODEL_SSNR_REF, "rlnSignalToNoiseRatio", "Spectral signal-to-noise ratio for a reference");
        EMDL::addLabel<double>(EMDL::MLMODEL_TAU2_FUDGE_FACTOR, "rlnTau2FudgeFactor", "Regularisation parameter with which estimates for the power in the references will be multiplied (T in original paper)");
        EMDL::addLabel<double>(EMDL::MLMODEL_TAU2_REF, "rlnReferenceTau2", "Spherical average of the estimated power in the signal of a reference");

        EMDL::addLabel<double>(EMDL::OPTIMISER_ACCURACY_ROT, "rlnOverallAccuracyRotations", "Overall accuracy of the rotational assignments (in degrees)");
        EMDL::addLabel<double>(EMDL::OPTIMISER_ACCURACY_TRANS, "rlnOverallAccuracyTranslations", "Overall accuracy of the translational assignments (in pixels)");
        EMDL::addLabel<double>(EMDL::OPTIMISER_ACCURACY_TRANS_ANGSTROM, "rlnOverallAccuracyTranslationsAngst", "Overall accuracy of the translational assignments (in Angstroms)");
        EMDL::addLabel<double>(EMDL::OPTIMISER_ADAPTIVE_FRACTION, "rlnAdaptiveOversampleFraction", "Fraction of the weights that will be oversampled in a second pass of the adaptive oversampling strategy");
        EMDL::addLabel<int>(EMDL::OPTIMISER_ADAPTIVE_OVERSAMPLING, "rlnAdaptiveOversampleOrder", "Order of the adaptive oversampling (0=no oversampling, 1= 2x oversampling; 2= 4x oversampling, etc)");
        EMDL::addLabel<int>(EMDL::OPTIMISER_AUTO_LOCAL_HP_ORDER, "rlnAutoLocalSearchesHealpixOrder", "Healpix order (before oversampling) from which autosampling procedure will use local angular searches");
        EMDL::addLabel<double>(EMDL::OPTIMISER_AVAILABLE_MEMORY, "rlnAvailableMemory", "Available memory per computing node (i.e. per MPI-process)");
        EMDL::addLabel<double>(EMDL::OPTIMISER_BEST_RESOL_THUS_FAR, "rlnBestResolutionThusFar", "The highest resolution that has been obtained in this optimization thus far");
        EMDL::addLabel<int>(EMDL::OPTIMISER_COARSE_SIZE, "rlnCoarseImageSize", "Current size of the images to be used in the first pass of the adaptive oversampling strategy (may be smaller than the original image size)");
        EMDL::addLabel<double>(EMDL::OPTIMISER_CHANGES_OPTIMAL_OFFSETS, "rlnChangesOptimalOffsets", "The average change in optimal translation in the last iteration (in pixels) ");
        EMDL::addLabel<double>(EMDL::OPTIMISER_CHANGES_OPTIMAL_ORIENTS, "rlnChangesOptimalOrientations", "The average change in optimal orientation in the last iteration (in degrees) ");
        EMDL::addLabel<double>(EMDL::OPTIMISER_CHANGES_OPTIMAL_CLASSES, "rlnChangesOptimalClasses", "The number of particles that changed their optimal clsas assignment in the last iteration");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_DATA_ARE_CTF_PHASE_FLIPPED, "rlnCtfDataArePhaseFlipped", "Flag to indicate that the input images have been phase-flipped");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED, "rlnCtfDataAreCtfPremultiplied", "Flag to indicate that the input images have been premultiplied with their CTF");
        EMDL::addLabel<std::string>(EMDL::OPTIMISER_DATA_STARFILE, "rlnExperimentalDataStarFile", "STAR file with metadata for the experimental images");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_DO_CORRECT_CTF, "rlnDoCorrectCtf", "Flag to indicate that CTF-correction should be performed");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_DO_CORRECT_MAGNIFICATION, "rlnDoCorrectMagnification", "Flag to indicate that (per-group) magnification correction should be performed");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_DO_CORRECT_NORM, "rlnDoCorrectNorm", "Flag to indicate that (per-image) normalisation-error correction should be performed");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_DO_CORRECT_SCALE, "rlnDoCorrectScale", "Flag to indicate that internal (per-group) intensity-scale correction should be performed");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_DO_EXTERNAL_RECONSTRUCT, "rlnDoExternalReconstruct", "Flag to indicate that the reconstruction will be performed outside relion_refine, e.g. for learned priors");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_DO_REALIGN_MOVIES, "rlnDoRealignMovies", "Flag to indicate that individual frames of movies are being re-aligned");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_DO_MAP, "rlnDoMapEstimation", "Flag to indicate that MAP estimation should be performed (otherwise ML estimation)");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_DO_SGD, "rlnDoStochasticGradientDescent", "Flag to indicate that SGD-optimisation should be performed (otherwise expectation maximisation)");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_DO_STOCHASTIC_EM, "rlnDoStochasticEM", "Flag to indicate that stochastic EM-optimisation should be performed (an alternative to SGD)");
        EMDL::addLabel<std::string>(EMDL::OPTIMISER_EXTERNAL_RECONS_DATA_REAL, "rlnExtReconsDataReal", "Name of the map with the real components of the input data array for the external reconstruction program");
        EMDL::addLabel<std::string>(EMDL::OPTIMISER_EXTERNAL_RECONS_DATA_IMAG, "rlnExtReconsDataImag", "Name of the map with the imaginary components of the input data array for the external reconstruction program");
        EMDL::addLabel<std::string>(EMDL::OPTIMISER_EXTERNAL_RECONS_WEIGHT, "rlnExtReconsWeight", "Name of the map with the input weight array for the external reconstruction program");
        EMDL::addLabel<std::string>(EMDL::OPTIMISER_EXTERNAL_RECONS_RESULT, "rlnExtReconsResult", "Name of the output reconstruction from the external reconstruction program");
        EMDL::addLabel<std::string>(EMDL::OPTIMISER_EXTERNAL_RECONS_NEWSTAR, "rlnExtReconsResultStarfile", "Name of the output STAR file with updated FSC or tau curves");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_FAST_SUBSETS, "rlnDoFastSubsetOptimisation", "Use subsets of the data in the earlier iterations to speed up convergence");
        EMDL::addLabel<int>(EMDL::OPTIMISER_SGD_INI_ITER, "rlnSgdInitialIterations", "Number of initial SGD iterations (at rlnSgdInitialResolution and with rlnSgdInitialSubsetSize)");
        EMDL::addLabel<int>(EMDL::OPTIMISER_SGD_FIN_ITER, "rlnSgdFinalIterations", "Number of final SGD iterations (at rlnSgdFinalResolution and with rlnSgdFinalSubsetSize)");
        EMDL::addLabel<int>(EMDL::OPTIMISER_SGD_INBETWEEN_ITER, "rlnSgdInBetweenIterations", "Number of SGD iteration in between the initial ones to the final ones (with linear interpolation of resolution and subset size)");
        EMDL::addLabel<double>(EMDL::OPTIMISER_SGD_INI_RESOL, "rlnSgdInitialResolution", "Resolution (in A) to use during the initial SGD iterations");
        EMDL::addLabel<double>(EMDL::OPTIMISER_SGD_FIN_RESOL, "rlnSgdFinalResolution", "Resolution (in A) to use during the final SGD iterations");
        EMDL::addLabel<int>(EMDL::OPTIMISER_SGD_INI_SUBSET_SIZE, "rlnSgdInitialSubsetSize", "Number of particles in a mini-batch (subset) during the initial SGD iterations");
        EMDL::addLabel<int>(EMDL::OPTIMISER_SGD_FIN_SUBSET_SIZE, "rlnSgdFinalSubsetSize", "Number of particles in a mini-batch (subset) during the final SGD iteration");
        EMDL::addLabel<double>(EMDL::OPTIMISER_SGD_MU, "rlnSgdMuFactor", "The mu-parameter that controls the momentum of the SGD gradients");
        EMDL::addLabel<double>(EMDL::OPTIMISER_SGD_SIGMA2FUDGE_INI, "rlnSgdSigma2FudgeInitial", "The variance of the noise will initially be multiplied with this value (larger than 1)");
        EMDL::addLabel<int>(EMDL::OPTIMISER_SGD_SIGMA2FUDGE_HALFLIFE, "rlnSgdSigma2FudgeHalflife", "After processing this many particles the multiplicative factor for the noise variance will have halved");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_SGD_SKIP_ANNNEAL, "rlnSgdSkipAnneal", "Option to switch off annealing of multiple references in SGD");
        EMDL::addLabel<int>(EMDL::OPTIMISER_SGD_SUBSET_SIZE, "rlnSgdSubsetSize", "The number of particles in the random subsets for SGD");
        EMDL::addLabel<int>(EMDL::OPTIMISER_SGD_WRITE_EVERY_SUBSET, "rlnSgdWriteEverySubset", "Every this many iterations the model is written to disk in SGD");
        EMDL::addLabel<int>(EMDL::OPTIMISER_SGD_MAX_SUBSETS, "rlnSgdMaxSubsets", "Stop SGD after doing this many subsets (possibly spanning more than 1 iteration)");
        EMDL::addLabel<double>(EMDL::OPTIMISER_SGD_STEPSIZE, "rlnSgdStepsize", "Stepsize in SGD updates)");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_DO_AUTO_REFINE, "rlnDoAutoRefine", "Flag to indicate that 3D auto-refine procedure is being used");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_DO_ONLY_FLIP_CTF_PHASES, "rlnDoOnlyFlipCtfPhases", "Flag to indicate that CTF-correction should only comprise phase-flipping");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_DO_SOLVENT_FLATTEN, "rlnDoSolventFlattening", "Flag to indicate that the references should be masked to set their solvent areas to a constant density");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_DO_SOLVENT_FSC, "rlnDoSolventFscCorrection", "Flag to indicate that the FSCs should be solvent-corrected during refinement");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_DO_SKIP_ALIGN, "rlnDoSkipAlign", "Flag to indicate that orientational (i.e. rotational and translational) searches will be omitted from the refinement, only marginalisation over classes will take place");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_DO_SKIP_ROTATE, "rlnDoSkipRotate", "Flag to indicate that rotational searches will be omitted from the refinement, only marginalisation over classes and translations will take place");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_DO_SPLIT_RANDOM_HALVES, "rlnDoSplitRandomHalves", "Flag to indicate that the data should be split into two completely separate, random halves");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_DO_ZERO_MASK, "rlnDoZeroMask", "Flag to indicate that the surrounding solvent area in the experimental particles will be masked to zeros (by default random noise will be used");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_FIX_SIGMA_NOISE, "rlnFixSigmaNoiseEstimates", "Flag to indicate that the estimates for the power spectra of the noise should be kept constant");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_FIX_SIGMA_OFFSET, "rlnFixSigmaOffsetEstimates", "Flag to indicate that the estimates for the stddev in the origin offsets should be kept constant");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_FIX_TAU, "rlnFixTauEstimates", "Flag to indicate that the estimates for the power spectra of the signal (i.e. the references) should be kept constant");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_HAS_CONVERGED, "rlnHasConverged", "Flag to indicate that the optimization has converged");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_HAS_HIGH_FSC_AT_LIMIT, "rlnHasHighFscAtResolLimit", "Flag to indicate that the FSC at the resolution limit is significant");
        EMDL::addLabel<int>(EMDL::OPTIMISER_HAS_LARGE_INCR_SIZE_ITER_AGO, "rlnHasLargeSizeIncreaseIterationsAgo", "How many iterations have passed since the last large increase in image size");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_DO_HELICAL_REFINE, "rlnDoHelicalRefine", "Flag to indicate that helical refinement should be performed");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_IGNORE_HELICAL_SYMMETRY, "rlnIgnoreHelicalSymmetry", "Flag to indicate that helical symmetry is ignored in 3D reconstruction");
        EMDL::addLabel<std::string>(EMDL::OPTIMISER_FOURIER_MASK, "rlnFourierMask", "Name of an FFTW-centred Fourier mask to be applied to the Projector for refinement.");
        EMDL::addLabel<double>(EMDL::OPTIMISER_HELICAL_TWIST_INITIAL, "rlnHelicalTwistInitial", "The intial helical twist (rotation per subunit) in degrees before refinement");
        EMDL::addLabel<double>(EMDL::OPTIMISER_HELICAL_RISE_INITIAL, "rlnHelicalRiseInitial", "The initial helical rise (translation per subunit) in Angstroms before refinement");
        EMDL::addLabel<double>(EMDL::OPTIMISER_HELICAL_Z_PERCENTAGE, "rlnHelicalCentralProportion", "Only expand this central fraction of the Z axis when imposing real-space helical symmetry");
        EMDL::addLabel<int>(EMDL::OPTIMISER_HELICAL_NSTART, "rlnNrHelicalNStart", "The N-number for an N-start helix");
        EMDL::addLabel<double>(EMDL::OPTIMISER_HELICAL_TUBE_INNER_DIAMETER, "rlnHelicalMaskTubeInnerDiameter", "Inner diameter of helical tubes in Angstroms (for masks of helical references and particles)");
        EMDL::addLabel<double>(EMDL::OPTIMISER_HELICAL_TUBE_OUTER_DIAMETER, "rlnHelicalMaskTubeOuterDiameter", "Outer diameter of helical tubes in Angstroms (for masks of helical references and particles)");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_HELICAL_SYMMETRY_LOCAL_REFINEMENT, "rlnHelicalSymmetryLocalRefinement", "Flag to indicate that local refinement of helical parameters should be performed");
        EMDL::addLabel<double>(EMDL::OPTIMISER_HELICAL_SIGMA_DISTANCE, "rlnHelicalSigmaDistance", "Sigma of distance along the helical tracks");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_HELICAL_KEEP_TILT_PRIOR_FIXED, "rlnHelicalKeepTiltPriorFixed", "Flag to indicate that helical tilt priors are kept fixed (at 90 degrees) in global angular searches");
        EMDL::addLabel<double>(EMDL::OPTIMISER_LOWRES_LIMIT_EXP, "rlnLowresLimitExpectation", "Low-resolution-limit (in Angstrom) for the expectation step");
        EMDL::addLabel<double>(EMDL::OPTIMISER_HIGHRES_LIMIT_EXP, "rlnHighresLimitExpectation", "High-resolution-limit (in Angstrom) for the expectation step");
        EMDL::addLabel<double>(EMDL::OPTIMISER_HIGHRES_LIMIT_SGD, "rlnHighresLimitSGD", "High-resolution-limit (in Angstrom) for Stochastic Gradient Descent");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_IGNORE_CTF_UNTIL_FIRST_PEAK, "rlnDoIgnoreCtfUntilFirstPeak", "Flag to indicate that the CTFs should be ignored until their first peak");
        EMDL::addLabel<int>(EMDL::OPTIMISER_INCR_SIZE, "rlnIncrementImageSize", "Number of Fourier shells to be included beyond the resolution where SSNR^MAP drops below 1");
        EMDL::addLabel<int>(EMDL::OPTIMISER_ITERATION_NO, "rlnCurrentIteration", "The number of the current iteration");
        EMDL::addLabel<std::string>(EMDL::OPTIMISER_LOCAL_SYMMETRY_FILENAME, "rlnLocalSymmetryFile", "Local symmetry description file containing list of masks and their operators");
        EMDL::addLabel<double>(EMDL::OPTIMISER_LOWRES_JOIN_RANDOM_HALVES, "rlnJoinHalvesUntilThisResolution", "Resolution (in Angstrom) to join the two random half-reconstructions to prevent their diverging orientations (for C-symmetries)");
        EMDL::addLabel<double>(EMDL::OPTIMISER_MAGNIFICATION_RANGE, "rlnMagnificationSearchRange", "Search range for magnification correction");
        EMDL::addLabel<double>(EMDL::OPTIMISER_MAGNIFICATION_STEP, "rlnMagnificationSearchStep", "Step size  for magnification correction");
        EMDL::addLabel<int>(EMDL::OPTIMISER_MAX_COARSE_SIZE, "rlnMaximumCoarseImageSize", "Maximum size of the images to be used in the first pass of the adaptive oversampling strategy (may be smaller than the original image size)");
        EMDL::addLabel<int>(EMDL::OPTIMISER_MAX_NR_POOL, "rlnMaxNumberOfPooledParticles", "Maximum number particles that are processed together to speed up calculations");
        EMDL::addLabel<std::string>(EMDL::OPTIMISER_MODEL_STARFILE, "rlnModelStarFile", "STAR file with metadata for the model that is being refined");
        EMDL::addLabel<std::string>(EMDL::OPTIMISER_MODEL_STARFILE2, "rlnModelStarFile2", "STAR file with metadata for the second model that is being refined (from random halves of the data)");
        EMDL::addLabel<int>(EMDL::OPTIMISER_NR_ITERATIONS, "rlnNumberOfIterations", "Maximum number of iterations to be performed");
        EMDL::addLabel<int>(EMDL::OPTIMISER_NR_ITER_WO_RESOL_GAIN, "rlnNumberOfIterWithoutResolutionGain", "Number of iterations that have passed without a gain in resolution");
        EMDL::addLabel<int>(EMDL::OPTIMISER_NR_ITER_WO_HIDDEN_VAR_CHANGES, "rlnNumberOfIterWithoutChangingAssignments", "Number of iterations that have passed without large changes in orientation and class assignments");
        EMDL::addLabel<std::string>(EMDL::OPTIMISER_OPTICS_STARFILE, "rlnOpticsStarFile", "STAR file with metadata for the optical groups (new as of version 3.1)");
        EMDL::addLabel<std::string>(EMDL::OPTIMISER_OUTPUT_ROOTNAME, "rlnOutputRootName", "Rootname for all output files (this may include a directory structure, which should then exist)");
        EMDL::addLabel<double>(EMDL::OPTIMISER_PARTICLE_DIAMETER, "rlnParticleDiameter", "Diameter of the circular mask to be applied to all experimental images (in Angstroms)");
        EMDL::addLabel<int>(EMDL::OPTIMISER_RADIUS_MASK_3D_MAP, "rlnRadiusMaskMap", "Radius of the spherical mask to be applied to all references (in Angstroms)");
        EMDL::addLabel<int>(EMDL::OPTIMISER_RADIUS_MASK_EXP_PARTICLES, "rlnRadiusMaskExpImages", "Radius of the circular mask to be applied to all experimental images (in Angstroms)");
        EMDL::addLabel<int>(EMDL::OPTIMISER_RANDOM_SEED, "rlnRandomSeed", "Seed (i.e. a number) for the random number generator");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_REFS_ARE_CTF_CORRECTED, "rlnRefsAreCtfCorrected", "Flag to indicate that the input references have been CTF-amplitude corrected");
        EMDL::addLabel<int>(EMDL::OPTIMISER_SMALLEST_CHANGES_OPT_CLASSES, "rlnSmallestChangesClasses", "Smallest changes thus far in the optimal class assignments (in numer of particles).");
        EMDL::addLabel<double>(EMDL::OPTIMISER_SMALLEST_CHANGES_OPT_OFFSETS, "rlnSmallestChangesOffsets", "Smallest changes thus far in the optimal offset assignments (in pixels).");
        EMDL::addLabel<double>(EMDL::OPTIMISER_SMALLEST_CHANGES_OPT_ORIENTS, "rlnSmallestChangesOrientations", "Smallest changes thus far in the optimal orientation assignments (in degrees).");
        EMDL::addLabel<std::string>(EMDL::OPTIMISER_SAMPLING_STARFILE, "rlnOrientSamplingStarFile", "STAR file with metadata for the orientational sampling");
        EMDL::addLabel<std::string>(EMDL::OPTIMISER_SOLVENT_MASK_NAME, "rlnSolventMaskName", "Name of an image that contains a (possibly soft) mask for the solvent area (values=0 for solvent, values =1 for protein)");
        EMDL::addLabel<std::string>(EMDL::OPTIMISER_SOLVENT_MASK2_NAME, "rlnSolventMask2Name", "Name of a secondary solvent mask (e.g. to flatten density inside an icosahedral virus)");
        EMDL::addLabel<std::string>(EMDL::OPTIMISER_TAU_SPECTRUM_NAME, "rlnTauSpectrumName", "Name of a STAR file that holds a tau2-spectrum");
        EMDL::addLabel<bool>(EMDL::OPTIMISER_USE_TOO_COARSE_SAMPLING, "rlnUseTooCoarseSampling", "Flag to indicate that the angular sampling on the sphere will be one step coarser than needed to speed up calculations");
        EMDL::addLabel<int>(EMDL::OPTIMISER_WIDTH_MASK_EDGE, "rlnWidthMaskEdge", "Width (in pixels) of the soft edge for spherical/circular masks to be used for solvent flattening");

        EMDL::addLabel<bool>(EMDL::ORIENT_FLIP, "rlnIsFlip", "Flag to indicate that an image should be mirrored");
        EMDL::addLabel<int>(EMDL::ORIENT_ID, "rlnOrientationsID", "ID (i.e. a unique number) for an orientation");
        EMDL::addLabel<double>(EMDL::ORIENT_ORIGIN_X, "rlnOriginX", "X-coordinate (in pixels) for the origin of rotation");
        EMDL::addLabel<double>(EMDL::ORIENT_ORIGIN_Y, "rlnOriginY", "Y-coordinate (in pixels) for the origin of rotation");
        EMDL::addLabel<double>(EMDL::ORIENT_ORIGIN_Z, "rlnOriginZ", "Z-coordinate (in pixels) for the origin of rotation");
        EMDL::addLabel<double>(EMDL::ORIENT_ORIGIN_X_PRIOR, "rlnOriginXPrior", "Center of the prior on the X-coordinate (in pixels) for the origin of rotation");
        EMDL::addLabel<double>(EMDL::ORIENT_ORIGIN_Y_PRIOR, "rlnOriginYPrior", "Center of the prior on the Y-coordinate (in pixels) for the origin of rotation");
        EMDL::addLabel<double>(EMDL::ORIENT_ORIGIN_Z_PRIOR, "rlnOriginZPrior", "Center of the prior on the Z-coordinate (in pixels) for the origin of rotation");

        EMDL::addLabel<double>(EMDL::ORIENT_ORIGIN_X_ANGSTROM, "rlnOriginXAngst", "X-coordinate (in Angstrom) for the origin of rotation");
        EMDL::addLabel<double>(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, "rlnOriginYAngst", "Y-coordinate (in Angstrom) for the origin of rotation");
        EMDL::addLabel<double>(EMDL::ORIENT_ORIGIN_Z_ANGSTROM, "rlnOriginZAngst", "Z-coordinate (in Angstrom) for the origin of rotation");
        EMDL::addLabel<double>(EMDL::ORIENT_ORIGIN_X_PRIOR_ANGSTROM, "rlnOriginXPriorAngst", "Center of the prior on the X-coordinate (in Angstrom) for the origin of rotation");
        EMDL::addLabel<double>(EMDL::ORIENT_ORIGIN_Y_PRIOR_ANGSTROM, "rlnOriginYPriorAngst", "Center of the prior on the Y-coordinate (in Angstrom) for the origin of rotation");
        EMDL::addLabel<double>(EMDL::ORIENT_ORIGIN_Z_PRIOR_ANGSTROM, "rlnOriginZPriorAngst", "Center of the prior on the Z-coordinate (in Angstrom) for the origin of rotation");

        EMDL::addLabel<double>(EMDL::ORIENT_ROT, "rlnAngleRot", "First Euler angle (rot, in degrees)");
        EMDL::addLabel<double>(EMDL::ORIENT_ROT_PRIOR, "rlnAngleRotPrior", "Center of the prior (in degrees) on the first Euler angle (rot)");
        EMDL::addLabel<double>(EMDL::ORIENT_ROT_PRIOR_FLIP_RATIO, "rlnAngleRotFlipRatio", "Flip ratio of bimodal rot prior (0~0.5, 0 means an ordinary prior, 0.5 means a perfect bimodal prior)");   // KThurber
        EMDL::addLabel<double>(EMDL::ORIENT_TILT, "rlnAngleTilt", "Second Euler angle (tilt, in degrees)");
        EMDL::addLabel<double>(EMDL::ORIENT_TILT_PRIOR, "rlnAngleTiltPrior", "Center of the prior (in degrees) on the second Euler angle (tilt)");
        EMDL::addLabel<double>(EMDL::ORIENT_PSI, "rlnAnglePsi", "Third Euler, or in-plane angle (psi, in degrees)");
        EMDL::addLabel<double>(EMDL::ORIENT_PSI_PRIOR, "rlnAnglePsiPrior", "Center of the prior (in degrees) on the third Euler angle (psi)");
        EMDL::addLabel<double>(EMDL::ORIENT_PSI_PRIOR_FLIP_RATIO, "rlnAnglePsiFlipRatio", "Flip ratio of bimodal psi prior (0~0.5, 0 means an ordinary prior, 0.5 means a perfect bimodal prior)");
        EMDL::addLabel<bool>(EMDL::ORIENT_PSI_PRIOR_FLIP, "rlnAnglePsiFlip", "Flag to indicate that psi prior angle has been flipped");  // KThurber

        EMDL::addLabel<double>(EMDL::PARTICLE_AUTOPICK_FOM, "rlnAutopickFigureOfMerit", "Autopicking FOM for a particle");
        EMDL::addLabel<int>(EMDL::PARTICLE_HELICAL_TUBE_ID, "rlnHelicalTubeID", "Helical tube ID for a helical segment");
        EMDL::addLabel<double>(EMDL::PARTICLE_HELICAL_TUBE_PITCH, "rlnHelicalTubePitch", "Cross-over distance for a helical segment (A)");
        EMDL::addLabel<double>(EMDL::PARTICLE_HELICAL_TRACK_LENGTH, "rlnHelicalTrackLength", "Distance (in pix) from the position of this helical segment to the starting point of the tube");
        EMDL::addLabel<double>(EMDL::PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM, "rlnHelicalTrackLengthAngst", "Distance (in A) from the position of this helical segment to the starting point of the tube");
        EMDL::addLabel<int>(EMDL::PARTICLE_CLASS, "rlnClassNumber", "Class number for which a particle has its highest probability");
        EMDL::addLabel<double>(EMDL::PARTICLE_DLL, "rlnLogLikeliContribution", "Contribution of a particle to the log-likelihood target function");
        EMDL::addLabel<int>(EMDL::PARTICLE_ID, "rlnParticleId", "ID (i.e. a unique number) for a particle");
        EMDL::addLabel<double>(EMDL::PARTICLE_FOM, "rlnParticleFigureOfMerit", "Developmental FOM for a particle");
        EMDL::addLabel<double>(EMDL::PARTICLE_KL_DIVERGENCE, "rlnKullbackLeiblerDivergence", "Kullback-Leibler divergence for a particle");
        EMDL::addAltLabel(EMDL::PARTICLE_KL_DIVERGENCE, "rlnKullbackLeibnerDivergence"); // wrong spelling for backwards compatibility
        EMDL::addLabel<int>(EMDL::PARTICLE_RANDOM_SUBSET, "rlnRandomSubset", "Random subset to which this particle belongs");
        EMDL::addLabel<int>(EMDL::PARTICLE_BEAM_TILT_CLASS, "rlnBeamTiltClass", "Beam-tilt class of a particle");
        EMDL::addLabel<std::string>(EMDL::PARTICLE_NAME, "rlnParticleName", "Name for a particle");
        EMDL::addLabel<std::string>(EMDL::PARTICLE_ORI_NAME, "rlnOriginalParticleName", "Original name for a particles");
        EMDL::addLabel<int>(EMDL::PARTICLE_NR_SIGNIFICANT_SAMPLES, "rlnNrOfSignificantSamples", "Number of orientational/class assignments (for a particle) with sign.probabilities in the 1st pass of adaptive oversampling"); /**< particle, Number of orientations contributing to weights*/
        EMDL::addLabel<int>(EMDL::PARTICLE_NR_FRAMES, "rlnNrOfFrames", "Number of movie frames that were collected for this particle");
        EMDL::addLabel<int>(EMDL::PARTICLE_NR_FRAMES_AVG, "rlnAverageNrOfFrames", "Number of movie frames that one averages over upon extraction of movie-particles");
        EMDL::addLabel<int>(EMDL::PARTICLE_MOVIE_RUNNING_AVG, "rlnMovieFramesRunningAverage", "Number of movie frames inside the running average that will be used for movie-refinement");
        EMDL::addLabel<double>(EMDL::PARTICLE_PMAX, "rlnMaxValueProbDistribution", "Maximum value of the (normalised) probability function for a particle"); /**< particle, Maximum value of probability distribution */
        EMDL::addLabel<int>(EMDL::PARTICLE_NUMBER, "rlnParticleNumber", "Number of particles");

        EMDL::addLabel<int>(EMDL::PIPELINE_JOB_COUNTER, "rlnPipeLineJobCounter", "Number of the last job in the pipeline");
        EMDL::addLabel<std::string>(EMDL::PIPELINE_NODE_NAME, "rlnPipeLineNodeName", "Name of a Node in the pipeline");
        EMDL::addLabel<int>(EMDL::PIPELINE_NODE_TYPE, "rlnPipeLineNodeType", "Type of a Node in the pipeline");
        EMDL::addLabel<std::string>(EMDL::PIPELINE_PROCESS_ALIAS, "rlnPipeLineProcessAlias", "Alias of a Process in the pipeline");
        EMDL::addLabel<std::string>(EMDL::PIPELINE_PROCESS_NAME, "rlnPipeLineProcessName", "Name of a Process in the pipeline");
        EMDL::addLabel<int>(EMDL::PIPELINE_PROCESS_TYPE, "rlnPipeLineProcessType", "Type of a Process in the pipeline");
        EMDL::addLabel<int>(EMDL::PIPELINE_PROCESS_STATUS, "rlnPipeLineProcessStatus", "Status of a Process in the pipeline (running, scheduled, finished or cancelled)");
        EMDL::addLabel<std::string>(EMDL::PIPELINE_EDGE_FROM, "rlnPipeLineEdgeFromNode", "Name of the origin of an edge");
        EMDL::addLabel<std::string>(EMDL::PIPELINE_EDGE_TO, "rlnPipeLineEdgeToNode", "Name of the to-Node in an edge");
        EMDL::addLabel<std::string>(EMDL::PIPELINE_EDGE_PROCESS, "rlnPipeLineEdgeProcess", "Name of the destination of an edge");

        EMDL::addLabel<double>(EMDL::POSTPROCESS_FINAL_RESOLUTION, "rlnFinalResolution", "Final estimated resolution after postprocessing (in Angstroms)");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_BFACTOR, "rlnBfactorUsedForSharpening", "Applied B-factor in the sharpening of the map");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_FRACTION_MOLWEIGHT, "rlnParticleBoxFractionMolecularWeight", "Fraction of protein voxels in the box, based on ordered molecular weight estimate, for calculating cisTEM-like part_FSC");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_FRACTION_SOLVENT_MASK, "rlnParticleBoxFractionSolventMask", "Fraction of protein voxels in the box, based on the solvent mask, for calculating cisTEM-like part_FSC");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_FSC_GENERAL, "rlnFourierShellCorrelation", "FSC value (of unspecified type, e.g. masked or unmasked)");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_FSC_TRUE, "rlnFourierShellCorrelationCorrected", "Final FSC value: i.e. after correction based on masking of randomized-phases maps");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_FSC_PART_MOLWEIGHT, "rlnFourierShellCorrelationParticleMolWeight", "CisTEM-like correction of unmasked FSCs, based on ordered molecular weight estimate");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_FSC_PART_FRACMASK, "rlnFourierShellCorrelationParticleMaskFraction", "CisTEM-like correction of unmasked FSCs, based on fraction of white pixels in solvent mask");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_FSC_MASKED, "rlnFourierShellCorrelationMaskedMaps", "FSC value after masking of the original maps");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_FSC_UNMASKED, "rlnFourierShellCorrelationUnmaskedMaps", "FSC value before masking of the original maps");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_FSC_RANDOM_MASKED, "rlnCorrectedFourierShellCorrelationPhaseRandomizedMaskedMaps", "FSC value after masking of the randomized-phases maps");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_AMPLCORR_MASKED, "rlnAmplitudeCorrelationMaskedMaps", "Correlation coefficient between amplitudes in Fourier shells of masked maps");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_AMPLCORR_UNMASKED, "rlnAmplitudeCorrelationUnmaskedMaps", "Correlation coefficient between amplitudes in Fourier shells of unmasked maps");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_DPR_MASKED, "rlnDifferentialPhaseResidualMaskedMaps", "Differential Phase Residual in Fourier shells of masked maps");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_DPR_UNMASKED, "rlnDifferentialPhaseResidualUnmaskedMaps", "Differential Phase Residual in Fourier shells of unmasked maps");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_GUINIER_FIT_INTERCEPT, "rlnFittedInterceptGuinierPlot", "The fitted intercept of the Guinier-plot");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_GUINIER_FIT_SLOPE, "rlnFittedSlopeGuinierPlot", "The fitted slope of the Guinier-plot");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_GUINIER_FIT_CORRELATION, "rlnCorrelationFitGuinierPlot", "The correlation coefficient of the fitted line through the Guinier-plot");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_GUINIER_VALUE_IN, "rlnLogAmplitudesOriginal", "Y-value for Guinier plot: the logarithm of the radially averaged amplitudes of the input map");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_GUINIER_VALUE_INVMTF, "rlnLogAmplitudesMTFCorrected", "Y-value for Guinier plot: the logarithm of the radially averaged amplitudes after MTF correction");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_GUINIER_VALUE_WEIGHTED, "rlnLogAmplitudesWeighted", "Y-value for Guinier plot: the logarithm of the radially averaged amplitudes after FSC-weighting");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_GUINIER_VALUE_SHARPENED, "rlnLogAmplitudesSharpened", "Y-value for Guinier plot: the logarithm of the radially averaged amplitudes after sharpening");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_GUINIER_VALUE_INTERCEPT, "rlnLogAmplitudesIntercept", "Y-value for Guinier plot: the fitted plateau of the logarithm of the radially averaged amplitudes");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_GUINIER_RESOL_SQUARED, "rlnResolutionSquared", "X-value for Guinier plot: squared resolution in 1/Angstrom^2");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_MOLWEIGHT, "rlnMolecularWeight", "Molecular weight of the ordered mass inside the box for calculating cisTEM-like part.FSC (in kDa)");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_MTF_VALUE, "rlnMtfValue", "Value of the detectors modulation transfer function (between 0 and 1)");
        EMDL::addLabel<double>(EMDL::POSTPROCESS_RANDOMISE_FROM, "rlnRandomiseFrom", "Resolution (in A) from which the phases are randomised in the postprocessing step");
        EMDL::addLabel<std::string>(EMDL::POSTPROCESS_UNFIL_HALFMAP1, "rlnUnfilteredMapHalf1", "Name of the unfiltered map from halfset 1");
        EMDL::addLabel<std::string>(EMDL::POSTPROCESS_UNFIL_HALFMAP2, "rlnUnfilteredMapHalf2", "Name of the unfiltered map from halfset 2");

        EMDL::addLabel<bool>(EMDL::SAMPLING_IS_3D, "rlnIs3DSampling", "Flag to indicate this concerns a 3D sampling ");
        EMDL::addLabel<bool>(EMDL::SAMPLING_IS_3D_TRANS, "rlnIs3DTranslationalSampling", "Flag to indicate this concerns a x,y,z-translational sampling ");
        EMDL::addLabel<int>(EMDL::SAMPLING_HEALPIX_ORDER, "rlnHealpixOrder", "Healpix order for the sampling of the first two Euler angles (rot, tilt) on the 3D sphere");
        EMDL::addLabel<int>(EMDL::SAMPLING_HEALPIX_ORDER_ORI, "rlnHealpixOrderOriginal", "Original healpix order for the sampling of the first two Euler angles (rot, tilt) on the 3D sphere");
        EMDL::addLabel<double>(EMDL::SAMPLING_LIMIT_TILT, "rlnTiltAngleLimit", "Values to which to limit the tilt angles (positive for keeping side views, negative for keeping top views)");
        EMDL::addLabel<double>(EMDL::SAMPLING_OFFSET_RANGE, "rlnOffsetRange", "Search range for the origin offsets (in Angstroms)");
        EMDL::addLabel<double>(EMDL::SAMPLING_OFFSET_STEP, "rlnOffsetStep", "Step size for the searches in the origin offsets (in Angstroms)");
        EMDL::addLabel<double>(EMDL::SAMPLING_OFFSET_RANGE_ORI, "rlnOffsetRangeOriginal", "Original search range for the origin offsets (in Angstroms)");
        EMDL::addLabel<double>(EMDL::SAMPLING_OFFSET_STEP_ORI, "rlnOffsetStepOriginal", "Original step size for the searches in the origin offsets (in Angstroms)");
        EMDL::addLabel<double>(EMDL::SAMPLING_HELICAL_OFFSET_STEP, "rlnHelicalOffsetStep", "Step size for the searches of offsets along helical axis (in Angstroms)");
        EMDL::addLabel<double>(EMDL::SAMPLING_PERTURB, "rlnSamplingPerturbInstance", "Random instance of the random perturbation on the orientational sampling");
        EMDL::addLabel<double>(EMDL::SAMPLING_PERTURBATION_FACTOR, "rlnSamplingPerturbFactor", "Factor for random perturbation on the orientational sampling (between 0 no perturbation and 1 very strong perturbation)");
        EMDL::addLabel<double>(EMDL::SAMPLING_PSI_STEP, "rlnPsiStep", "Step size (in degrees) for the sampling of the in-plane rotation angle (psi)");
        EMDL::addLabel<double>(EMDL::SAMPLING_PSI_STEP_ORI, "rlnPsiStepOriginal", "Original step size (in degrees) for the sampling of the in-plane rotation angle (psi)");
        EMDL::addLabel<std::string>(EMDL::SAMPLING_SYMMETRY, "rlnSymmetryGroup", "Symmetry group (e.g., C1, D7, I2, I5, etc.)");

        EMDL::addLabel<int>(EMDL::SCHEDULE_EDGE_NUMBER, "rlnScheduleEdgeNumber", "Numbered index of an edge inside a Schedule");
        EMDL::addLabel<std::string>(EMDL::SCHEDULE_EDGE_INPUT, "rlnScheduleEdgeInputNodeName" , "Name of the input Node for a schedule Edge");
        EMDL::addLabel<std::string>(EMDL::SCHEDULE_EDGE_OUTPUT, "rlnScheduleEdgeOutputNodeName", "Name of the output Node for a schedule Edge");
        EMDL::addLabel<bool>(EMDL::SCHEDULE_EDGE_IS_FORK, "rlnScheduleEdgeIsFork", "Flag to indicate that this Edge is a Fork, dependent on a Boolean Schedule variable");
        EMDL::addLabel<std::string>(EMDL::SCHEDULE_EDGE_OUTPUT_TRUE, "rlnScheduleEdgeOutputNodeNameIfTrue", "Name of the output Node for a schedule Fork if the associated Boolean is True");
        EMDL::addLabel<std::string>(EMDL::SCHEDULE_EDGE_BOOLEAN, "rlnScheduleEdgeBooleanVariable", "Name of the associated Boolean variable if this Edge is a Fork");
        EMDL::addLabel<std::string>(EMDL::SCHEDULE_GENERAL_CURRENT_NODE, "rlnScheduleCurrentNodeName", "Name of the current Node for this Schedule");
        EMDL::addLabel<std::string>(EMDL::SCHEDULE_GENERAL_ORIGINAL_START_NODE, "rlnScheduleOriginalStartNodeName", "Name of the original starting Node for this Schedule");
        EMDL::addLabel<std::string>(EMDL::SCHEDULE_GENERAL_EMAIL, "rlnScheduleEmailAddress", "Email address to send message when Schedule finishes");
        EMDL::addLabel<std::string>(EMDL::SCHEDULE_GENERAL_NAME, "rlnScheduleName", "Name for this Schedule");
        EMDL::addLabel<std::string>(EMDL::SCHEDULE_JOB_NAME, "rlnScheduleJobName", "Name of a Job in a Schedule");
        EMDL::addLabel<std::string>(EMDL::SCHEDULE_JOB_ORI_NAME, "rlnScheduleJobNameOriginal", "Original name of a Job in a Schedule");
        EMDL::addLabel<std::string>(EMDL::SCHEDULE_JOB_MODE, "rlnScheduleJobMode", "Mode on how to execute a Job");
        EMDL::addLabel<bool>(EMDL::SCHEDULE_JOB_HAS_STARTED, "rlnScheduleJobHasStarted", "Flag to indicate whether a Job has started already in the execution of the Schedule");
        EMDL::addLabel<std::string>(EMDL::SCHEDULE_OPERATOR_NAME,   "rlnScheduleOperatorName", "Name of a Boolean operator in the Schedule");
        EMDL::addLabel<std::string>(EMDL::SCHEDULE_OPERATOR_TYPE,   "rlnScheduleOperatorType", "Type of an operator in the Schedule");
        EMDL::addLabel<std::string>(EMDL::SCHEDULE_OPERATOR_INPUT1, "rlnScheduleOperatorInput1", "Name of the 1st input to the operator");
        EMDL::addLabel<std::string>(EMDL::SCHEDULE_OPERATOR_INPUT2, "rlnScheduleOperatorInput2", "Name of the 2nd input to the operator");
        EMDL::addLabel<std::string>(EMDL::SCHEDULE_OPERATOR_OUTPUT, "rlnScheduleOperatorOutput", "Name of the output variable on which this operator acts");
        EMDL::addLabel<std::string>(EMDL::SCHEDULE_VAR_BOOL_NAME,   "rlnScheduleBooleanVariableName", "Name of a Boolean variable in the Schedule");
        EMDL::addLabel<bool>(EMDL::SCHEDULE_VAR_BOOL_VALUE, "rlnScheduleBooleanVariableValue", "Value of a Boolean variable in the Schedule");
        EMDL::addLabel<bool>(EMDL::SCHEDULE_VAR_BOOL_ORI_VALUE, "rlnScheduleBooleanVariableResetValue", "Value which a Boolean variable will take upon a reset");
        EMDL::addLabel<std::string>(EMDL::SCHEDULE_VAR_FLOAT_NAME, "rlnScheduleFloatVariableName", "Name of a Float variable in the Schedule");
        EMDL::addLabel<double>(EMDL::SCHEDULE_VAR_FLOAT_VALUE, "rlnScheduleFloatVariableValue", "Value of a Float variable in the Schedule");
        EMDL::addLabel<double>(EMDL::SCHEDULE_VAR_FLOAT_ORI_VALUE, "rlnScheduleFloatVariableResetValue", "Value which a Float variable will take upon a reset");
        EMDL::addLabel<std::string>(EMDL::SCHEDULE_VAR_STRING_NAME, "rlnScheduleStringVariableName", "Name of a String variable in the Schedule");
        EMDL::addLabel<std::string>(EMDL::SCHEDULE_VAR_STRING_VALUE, "rlnScheduleStringVariableValue", "Value of a String variable in the Schedule");
        EMDL::addLabel<std::string>(EMDL::SCHEDULE_VAR_STRING_ORI_VALUE, "rlnScheduleStringVariableResetValue", "Value which a String variable will take upon a reset");

        EMDL::addLabel<int>(EMDL::SELECTED, "rlnSelected", "Flag whether an entry in a metadatatable is selected (1) in the viewer or not (0)");
        EMDL::addLabel<double>(EMDL::SELECT_PARTICLES_ZSCORE, "rlnParticleSelectZScore", "Sum of Z-scores from particle_select. High Z-scores are likely to be outliers.");
        EMDL::addLabel<int>(EMDL::SORTED_IDX, "rlnSortedIndex", "Index of a metadata entry after sorting (first sorted index is 0).");
        EMDL::addLabel<std::string>(EMDL::STARFILE_MOVIE_PARTICLES, "rlnStarFileMovieParticles", "Filename of a STAR file with movie-particles in it");
        EMDL::addLabel<double>(EMDL::PERFRAME_CUMULATIVE_WEIGHT, "rlnPerFrameCumulativeWeight", "Sum of the resolution-dependent relative weights from the first frame until the given frame");
        EMDL::addLabel<double>(EMDL::PERFRAME_RELATIVE_WEIGHT, "rlnPerFrameRelativeWeight", "The resolution-dependent relative weights for a given frame");

        EMDL::addLabel<double>(EMDL::RESOLUTION, "rlnResolution", "Resolution (in 1/Angstroms)");
        EMDL::addLabel<double>(EMDL::RESOLUTION_ANGSTROM, "rlnAngstromResolution", "Resolution (in Angstroms)");
        EMDL::addLabel<double>(EMDL::RESOLUTION_INVPIXEL, "rlnResolutionInversePixel", "Resolution (in 1/pixel, Nyquist = 0.5)");
        EMDL::addLabel<int>(EMDL::SPECTRAL_IDX, "rlnSpectralIndex", "Spectral index (i.e. distance in pixels to the origin in Fourier space) ");

        EMDL::addLabel<void>(EMDL::UNKNOWN_LABEL, "rlnUnknownLabel", "NON-RELION label: values will be ignored, yet maintained in the STAR file.");
    }

    ~StaticInitialization() {}

    friend class EMDL;

};

#endif
