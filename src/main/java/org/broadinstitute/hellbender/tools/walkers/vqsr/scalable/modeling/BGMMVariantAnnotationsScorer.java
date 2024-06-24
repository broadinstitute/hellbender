package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling;

import org.apache.commons.lang3.NotImplementedException;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5LibException;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.utils.HDF5Utils;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.LabeledVariantAnnotationsData;
import org.broadinstitute.hellbender.utils.clustering.BayesianGaussianMixtureModeller;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.Serializable;
import java.util.List;

// TODO this is just a stub, will be fleshed out in a separate PR
public final class BGMMVariantAnnotationsScorer implements VariantAnnotationsScorer, Serializable {

    private static final long serialVersionUID = 1L;

    public static final String BGMM_SCORER_SER_SUFFIX = ".bgmmScorer.ser";

    public BGMMVariantAnnotationsScorer(final List<String> annotationNames,
                                        final BGMMVariantAnnotationsModel.Preprocesser preprocesser,
                                        final BayesianGaussianMixtureModeller bgmm) {
        throw new NotImplementedException("BGMM module will be implemented in separate PR.");
    }

    @Override
    public void score(final File inputAnnotationsFile,
                      final File outputScoresFile) {
        throw new NotImplementedException("BGMM module will be implemented in separate PR.");
    }

    public double[][] preprocess(final double[][] annotations) {
        throw new NotImplementedException("BGMM module will be implemented in separate PR.");
    }

    public void serialize(final File scorerFile) {
        throw new NotImplementedException("BGMM module will be implemented in separate PR.");
    }

    public static BGMMVariantAnnotationsScorer deserialize(final File scorerFile) {
        throw new NotImplementedException("BGMM module will be implemented in separate PR.");
    }

    // TODO clean this up, copy more fields
    public static void preprocessAnnotationsWithBGMMAndWriteHDF5(final List<String> annotationNames,
                                                                 final String outputPrefix,
                                                                 final File labeledTrainingAndVariantTypeAnnotationsFile,
                                                                 final Logger logger) {
        final double[][] rawAnnotations = LabeledVariantAnnotationsData.readAnnotations(labeledTrainingAndVariantTypeAnnotationsFile);
        final BGMMVariantAnnotationsScorer scorer = BGMMVariantAnnotationsScorer.deserialize(new File(outputPrefix + BGMM_SCORER_SER_SUFFIX));
        final double[][] preprocessedAnnotations = scorer.preprocess(rawAnnotations);
        final File outputPreprocessedAnnotationsFile = new File(outputPrefix + ".annot.pre.hdf5");
        try (final HDF5File hdf5File = new HDF5File(outputPreprocessedAnnotationsFile, HDF5File.OpenMode.CREATE)) {
            IOUtils.canReadFile(hdf5File.getFile());
            hdf5File.makeStringArray("/data/annotation_names", annotationNames.toArray(new String[0]));
            HDF5Utils.writeChunkedDoubleMatrix(hdf5File, "/data/annotations", preprocessedAnnotations, HDF5Utils.MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX / 16);
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during writing of preprocessed annotations (%s). Output file at %s may be in a bad state.",
                    exception, outputPreprocessedAnnotationsFile.getAbsolutePath()));
        }
        logger.info(String.format("Preprocessed annotations written to %s.", outputPreprocessedAnnotationsFile.getAbsolutePath()));
    }
}
