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

    class StaticInitialization;  // A class for static initialization

    static std::map<const EMDLabel,    const LabelData>   data;
    static std::map<const std::string, const EMDLabel>    labels;
    static std::map<const std::string, const std::string> definitions;
    static StaticInitialization initialization;

    template <typename T>
    static void addLabel(const EMDLabel label, const std::string &name, const std::string &definition = "undocumented");
    static void addAltLabel(const EMDLabel label, const std::string &name);

    friend class StaticInitialization;

};

#endif
