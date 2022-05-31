"""Helper functions for test modules."""
from pydicom.dataset import Dataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import generate_uid
import numpy


# Constants for helper methods below
SLICE_Z = [-20, -10, 0, 10, 20]
NUM_SLICES = len(SLICE_Z)
STUDY_iUID = generate_uid()
CT_SERIES_iUID = generate_uid()
CT_iUID = [generate_uid() for i in range(NUM_SLICES)]
DOSE_SERIES_iUID = generate_uid()
PLAN_SERIES_iUID = generate_uid()
PLAN_iUID = generate_uid()
SS_SERIES_iUID = generate_uid()
SS_iUID = generate_uid()
FoR_UID = generate_uid()
DOSE_iUID = generate_uid()

# Position, orientation
CT_ipp = [2, 12, -20]
CT_iop = [0, -1, 0, 1, 0, 0]  # head-first decubitus left
DOSE_ipp = CT_ipp
DOSE_iop = CT_iop

RTStructureSetStorage = "1.2.840.10008.5.1.4.1.1.481.3"
CTImageStorage = "1.2.840.10008.5.1.4.1.1.2"
RTDoseStorage = "1.2.840.10008.5.1.4.1.1.481.2"


# Helper methods for testing DVH calc, particularly decubitus
def _basic_file_meta(class_UID):
    """Create basic file meta for DICOM dataset."""
    # File meta info data elements
    file_meta = FileMetaDataset()
    file_meta.FileMetaInformationVersion = b'\x00\x01'
    file_meta.MediaStorageSOPClassUID = class_UID
    file_meta.MediaStorageSOPInstanceUID = generate_uid()
    file_meta.TransferSyntaxUID = '1.2.840.10008.1.2'
    return file_meta


def fake_rtdose():
    """Create a fake RT Dose DICOM dataset."""
    # Main data elements
    file_meta = _basic_file_meta(RTDoseStorage)
    ds = Dataset()
    ds.SOPClassUID = RTDoseStorage
    ds.SOPInstanceUID = DOSE_iUID
    ds.StudyDate = '20220101'
    ds.Modality = 'RTDOSE'
    ds.PatientName = 'Decubitus'
    ds.PatientID = 'Decubitus'
    ds.PatientBirthDate = '18000101'
    ds.PatientSex = 'O'

    ds.SliceThickness = None
    ds.StudyInstanceUID = STUDY_iUID
    ds.SeriesInstanceUID = DOSE_SERIES_iUID
    ds.StudyID = 'none'
    ds.SeriesNumber = 1
    ds.ImagePositionPatient = DOSE_ipp
    ds.ImageOrientationPatient = DOSE_iop
    ds.FrameOfReferenceUID = FoR_UID
    ds.PositionReferenceIndicator = ''

    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = 'MONOCHROME2'
    ds.NumberOfFrames = 5
    ds.FrameIncrementPointer = (0x3004, 0x000c)  # GridFrameOffsetVector
    ds.Rows = 7
    ds.Columns = 8
    ds.PixelSpacing = [1.0, 2.0]  # (between rows, between cols)
    ds.BitsAllocated = 32
    ds.BitsStored = 32
    ds.HighBit = 31
    ds.PixelRepresentation = 0
    ds.DoseUnits = 'GY'
    ds.DoseType = 'PHYSICAL'
    ds.DoseSummationType = 'PLAN'
    ds.GridFrameOffsetVector = [0.0, 10.0, 20.0, 30.0, 40.0]

    # Referenced RT Plan Sequence
    refd_rt_plan_sequence = Sequence()
    ds.ReferencedRTPlanSequence = refd_rt_plan_sequence

    # Referenced RT Plan Sequence: Referenced RT Plan 1
    refd_rt_plan1 = Dataset()
    refd_rt_plan1.ReferencedSOPClassUID = '1.2.840.10008.5.1.4.1.1.481.5'
    refd_rt_plan1.ReferencedSOPInstanceUID = PLAN_iUID
    refd_rt_plan_sequence.append(refd_rt_plan1)

    arr = numpy.zeros(
        (ds.NumberOfFrames, ds.Rows, ds.Columns), dtype=numpy.uint32
    )

    ds.DoseGridScaling = 0.00001  # take to near integer cGy in range of < 100
    # Add small amount so values stay above their int value when dose scaled
    rounding_guard = 1
    for zindex in range(len(arr)):
        for row in range(len(arr[0])):
            for col in range(len(arr[0][0])):
                arr[zindex, row, col] = 1000 * (
                    zindex*10 + (row > 2)*row + (col > 2) * (col)
                ) + rounding_guard

    # Three middle slices after above math:
    # # Shown with location of contours:
    # (3, 14.5), (7, 14.5), (7, 17.5), (3, 17.5)  (2 voxels across x, 3 down y)
    #     2 mm x 1 mm pixel, across (cols, X) center of pixels from 2 to 16
    #                         down (rows, Y) centers from 12 to 18
    # volume = (2*2) x (3*1) x (3*10) = 360 mm^3 = 0.36 cm^3
    #       X=2   4   6   8  10  12  14  16
    # Y=12  [10, 10, 10, 13, 14, 15, 16, 17]
    #       [10, 10, 10, 13, 14, 15, 16, 17]
    #       [10, 10, 10, 13, 14, 15, 16, 17]
    #           |-------|
    #   15  [13, 13, 13, 16, 17, 18, 19, 20]
    #       [14, 14, 14, 17, 18, 19, 20, 21]
    #       [15, 15, 15, 18, 19, 20, 21, 22]
    #           |-------|
    # Y=18  [16, 16, 16, 19, 20, 21, 22, 23]

    #       X=2   4   6   8  10  12  14  16
    # Y=12  [20, 20, 20, 23, 24, 25, 26, 27]
    #       [20, 20, 20, 23, 24, 25, 26, 27]
    #       [20, 20, 20, 23, 24, 25, 26, 27]
    #           |-------|
    #       [23, 23, 23, 26, 27, 28, 29, 30]
    #       [24, 24, 24, 27, 28, 29, 30, 31]
    #       [25, 25, 25, 28, 29, 30, 31, 32]
    #           |-------|
    # Y=18  [26, 26, 26, 29, 30, 31, 32, 33]

    #       X=2   4   6   8  10  12  14  16
    # Y=12  [30, 30, 30, 33, 34, 35, 36, 37]
    #       [30, 30, 30, 33, 34, 35, 36, 37]
    #       [30, 30, 30, 33, 34, 35, 36, 37]
    #           |-------|
    #       [33, 33, 33, 36, 37, 38, 39, 40]
    #       [34, 34, 34, 37, 38, 39, 40, 41]
    #       [35, 35, 35, 38, 39, 40, 41, 42]
    #           |-------|
    # Y=18  [36, 36, 36, 39, 40, 41, 42, 43]

    ds.PixelData = arr.tobytes()
    ds.file_meta = file_meta
    ds.is_implicit_VR = True
    ds.is_little_endian = True

    return ds


def fake_ss():
    """Create a fake RT Structure Set DICOM dataset."""
    file_meta = _basic_file_meta(RTStructureSetStorage)
    # Main data elements
    ds = Dataset()
    ds.SpecificCharacterSet = 'ISO_IR 192'
    ds.InstanceCreationDate = '20220101'
    ds.SOPClassUID = RTStructureSetStorage
    ds.SOPInstanceUID = SS_iUID
    ds.StudyDate = '20220101'
    ds.Modality = 'RTSTRUCT'

    ds.PatientName = 'Decubitus'
    ds.PatientID = 'Decubitus'
    ds.PatientBirthDate = '18000101'
    ds.PatientSex = 'O'
    ds.StudyInstanceUID = STUDY_iUID
    ds.SeriesInstanceUID = SS_SERIES_iUID
    ds.StudyID = '1'
    ds.SeriesNumber = '10'
    ds.StructureSetLabel = 'CT_1'
    ds.StructureSetDate = '20220101'

    # Referenced Frame of Reference Sequence
    refd_frame_of_ref_sequence = Sequence()
    ds.ReferencedFrameOfReferenceSequence = refd_frame_of_ref_sequence

    # Referenced Frame of Reference Sequence: Referenced Frame of Reference 1
    refd_frame_of_ref1 = Dataset()
    refd_frame_of_ref1.FrameOfReferenceUID = FoR_UID

    # RT Referenced Study Sequence
    rt_refd_study_sequence = Sequence()
    refd_frame_of_ref1.RTReferencedStudySequence = rt_refd_study_sequence

    # RT Referenced Study Sequence: RT Referenced Study 1
    rt_refd_study1 = Dataset()
    rt_refd_study1.ReferencedSOPClassUID = '1.2.840.10008.5.1.4.1.1.481.3'
    rt_refd_study1.ReferencedSOPInstanceUID = STUDY_iUID

    # RT Referenced Series Sequence
    rt_refd_series_sequence = Sequence()
    rt_refd_study1.RTReferencedSeriesSequence = rt_refd_series_sequence

    # RT Referenced Series Sequence: RT Referenced Series 1
    rt_refd_series1 = Dataset()
    rt_refd_series1.SeriesInstanceUID = CT_SERIES_iUID

    # Contour Image Sequence
    contour_image_sequence = Sequence()
    rt_refd_series1.ContourImageSequence = contour_image_sequence

    # Contour Image Sequence: Contour Image 1
    contour_image1 = Dataset()
    contour_image1.ReferencedSOPClassUID = CTImageStorage
    contour_image1.ReferencedSOPInstanceUID = CT_iUID
    contour_image_sequence.append(contour_image1)

    rt_refd_series_sequence.append(rt_refd_series1)
    rt_refd_study_sequence.append(rt_refd_study1)
    refd_frame_of_ref_sequence.append(refd_frame_of_ref1)

    # Structure Set ROI Sequence
    structure_set_roi_sequence = Sequence()
    ds.StructureSetROISequence = structure_set_roi_sequence

    # Structure Set ROI Sequence: Structure Set ROI 1
    structure_set_roi1 = Dataset()
    structure_set_roi1.ROINumber = 1
    structure_set_roi1.ReferencedFrameOfReferenceUID = FoR_UID
    structure_set_roi1.ROIName = 'PTV'
    structure_set_roi1.ROIGenerationAlgorithm = 'MANUAL'
    structure_set_roi_sequence.append(structure_set_roi1)

    # ROI Contour Sequence
    #    Contour Sequence
    #        Contour Image Sequence (optional, references the CT slices)
    #        Contour Geometric Type (e.g. CLOSED_PLANAR)
    #        Number of Contour Points
    #        Contour Data
    #    Referenced ROI Number

    roi_contour_sequence = Sequence()
    ds.ROIContourSequence = roi_contour_sequence

    # ROI Contour Sequence: ROI Contour 1
    roi_contour1 = Dataset()
    roi_contour1.ROIDisplayColor = [255, 128, 0]

    # Contour Sequence
    contour_sequence = Sequence()
    roi_contour1.ContourSequence = contour_sequence

    for i in range(1, NUM_SLICES-1):  # contours on all but first, last slice
        contour = Dataset()
        contour.ContourGeometricType = 'CLOSED_PLANAR'

        z = SLICE_Z[i]
        contour.ContourData = [3, 14.5, z, 7, 14.5, z, 7, 17.5, z, 3, 17.5, z]
        contour.NumberOfContourPoints = len(contour.ContourData) / 3
        contour_sequence.append(contour)

    roi_contour1.ReferencedROINumber = 1
    roi_contour_sequence.append(roi_contour1)

    ds.file_meta = file_meta
    ds.is_implicit_VR = True
    ds.is_little_endian = True

    return ds
