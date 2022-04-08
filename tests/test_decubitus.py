import unittest
from pydicom.dataset import Dataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import generate_uid, PYDICOM_IMPLEMENTATION_UID
import numpy
from dicompylercore.dvhcalc import get_dvh


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


def basic_file_meta(class_UID):
    # File meta info data elements
    file_meta = FileMetaDataset()
    file_meta.FileMetaInformationVersion = b'\x00\x01'
    file_meta.MediaStorageSOPClassUID = class_UID
    file_meta.MediaStorageSOPInstanceUID = generate_uid()
    file_meta.TransferSyntaxUID = '1.2.840.10008.1.2'
    file_meta.ImplementationClassUID = PYDICOM_IMPLEMENTATION_UID
    return file_meta


def fake_rtdose():
    # Main data elements
    file_meta = basic_file_meta(RTDoseStorage)
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
    file_meta = basic_file_meta(RTStructureSetStorage)
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


class TestDVHCalcDecubitus(unittest.TestCase):
    """Unit tests for DVH calculation in decubitus orientations."""

    def setUp(self):
        self.ss = fake_ss()
        self.dose = fake_rtdose()

    def test_nondecub(self):
        """Test that DVH is calculated correctly for standard orientation."""
        self.dose.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
        dvh = get_dvh(self.ss, self.dose, 1)
        diffl = dvh.differential
        # Counts are normalized to total, and to volume,
        # So undo that here for test dose grid.
        # 18=num dose voxels inside struct; 0.36=volume
        got_counts = diffl.counts * 18 / 0.36
        expected_counts = [0]*13 + [2, 2, 2, 0, 0, 0, 0, 0, 0, 0,
                                    2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2]
        assert numpy.all(numpy.isclose(got_counts, expected_counts))

    def test_HF_decubitus_left(self):
        """Test DVH for head-first decubitus left orientation."""
        # Keep same dose grid as std orientation but pixel-spacing in X, Y same
        # For this case, use iop=[0, -1, 0, 1, 0, 0] Head first decubitus left
        # Then X = r * dr + ipp[0]
        #  and Y = -c * dc + ipp[1]
        # (https://nipy.org/nibabel/dicom/dicom_orientation.html
        # #dicom-affine-formula)
        # Change ipp y of y to new max of 19 for similar y range
        # Below show contours box of (3, 14.5) - (7, 17.5) on dose grid
        #       Y=19 18  17                  12
        # X=2   [10, 10, 10, 13, 14, 15, 16, 17],
        #               |-----------|
        #   4   [10, 10, 10, 13, 14, 15, 16, 17]
        #   6   [10, 10, 10, 13, 14, 15, 16, 17]
        #               |-----------|
        #   8   [13, 13, 13, 16, 17, 18, 19, 20]
        #  10   [14, 14, 14, 17, 18, 19, 20, 21]
        #  12   [15, 15, 15, 18, 19, 20, 21, 22]
        #  14   [16, 16, 16, 19, 20, 21, 22, 23]]

        #       Y=19 18  17                  12
        # X=2   [20, 20, 20, 23, 24, 25, 26, 27]
        #               |-----------|
        #   4   [20, 20, 20, 23, 24, 25, 26, 27]
        #   6   [20, 20, 20, 23, 24, 25, 26, 27]
        #               |-----------|
        #   8   [23, 23, 23, 26, 27, 28, 29, 30]
        #  10   [24, 24, 24, 27, 28, 29, 30, 31]
        #  12   [25, 25, 25, 28, 29, 30, 31, 32]
        #  14   [...]

        #       Y=19 18  17                  12
        # X=2   [30, 30, 30, 33, 34, 35, 36, 37]
        #               |-----------|
        #   4   [30, 30, 30, 33, 34, 35, 36, 37]
        #   6   [30, 30, 30, 33, 34, 35, 36, 37]
        #               |-----------|
        #   8   [33, 33, 33, 36, 37, 38, 39, 40]
        #   10  [34, 34, 34, 37, 38, 39, 40, 41]
        #   12  [35, 35, 35, 38, 39, 40, 41, 42]
        # X=14  [36, 36, 36, 39, 40, 41, 42, 43]

        #                          10       13 14                20
        expected_counts = [0]*10 + [2, 0, 0, 2, 2, 0, 0, 0, 0, 0, 2, 0, 0,
                                    2, 2, 0, 0, 0, 0, 0, 2, 0, 0, 2, 2]
        #                          23 24                30       33 34
        self.dose.ImagePositionPatient = [2, 19, -20]  # X Y Z top left
        self.dose.PixelSpacing = [2.0, 1.0]  # between Rows, Columns
        dvh = get_dvh(self.ss, self.dose, 1)
        diffl = dvh.differential
        # Counts are normalized to total, and to volume,
        # So undo that here for test dose grid.
        # 18=num dose voxels inside struct; 0.36=volume
        got_counts = diffl.counts * 18 / 0.36
        assert numpy.all(numpy.isclose(got_counts, expected_counts))

    def test_HF_decubitus_left_structure_extents(self):
        """Test DVH for HF decubitus Lt orientation structure_extents used."""
        # Repeat test_HF_decubitus_left but with use_structure_extents
        #                          10       13 14                20
        expected_counts = [0]*10 + [2, 0, 0, 2, 2, 0, 0, 0, 0, 0, 2, 0, 0,
                                    2, 2, 0, 0, 0, 0, 0, 2, 0, 0, 2, 2]
        #                          23 24                30       33 34
        self.dose.ImagePositionPatient = [2, 19, -20]  # X Y Z top left
        self.dose.PixelSpacing = [2.0, 1.0]  # between Rows, Columns
        dvh = get_dvh(self.ss, self.dose, 1, use_structure_extents=True)
        diffl = dvh.differential
        # Counts are normalized to total, and to volume,
        # So undo that here for test dose grid.
        # 18=num dose voxels inside struct; 0.36=volume
        got_counts = diffl.counts * 18 / 0.36
        assert numpy.all(numpy.isclose(got_counts, expected_counts))

    def test_HF_decubitus_right(self):
        """Test DVH for head-first decubitus right orientation."""
        # Keep same dose grid as std orientation

        self.dose.ImageOrientationPatient = [0, 1, 0, -1, 0, 0]
        self.dose.PixelSpacing = [2.0, 1.0]  # between Rows, Columns
        # original ipp = [2, 12, -20]
        # Then X = -r * dr + ipp[0], X decreases down the rows
        #  and Y = c * dc + ipp[1], Y increases across cols
        # (https://nipy.org/nibabel/dicom/dicom_orientation.html
        # #dicom-affine-formula)
        # Change ipp y of X to new max of 14 for similar y range
        self.dose.ImagePositionPatient = [14, 12, -20]  # X Y Z top left
        # Below show contours box of (3, 14.5) - (7, 17.5) on dose grid
        #       Y=12 13  14  15  16  17  18  19
        # X=14  [10, 10, 10, 13, 14, 15, 16, 17],
        #   12  [10, 10, 10, 13, 14, 15, 16, 17]
        #   10  [10, 10, 10, 13, 14, 15, 16, 17]
        #    8  [13, 13, 13, 16, 17, 18, 19, 20]
        #                   | ----------|
        #    6  [14, 14, 14, 17, 18, 19, 20, 21]
        #    4  [15, 15, 15, 18, 19, 20, 21, 22]
        #                   | ----------|
        #    2  [16, 16, 16, 19, 20, 21, 22, 23]]

        #       Y=12 13  14                  19
        # X=14  [20, 20, 20, 23, 24, 25, 26, 27]
        #   12  [20, 20, 20, 23, 24, 25, 26, 27]
        #   10  [20, 20, 20, 23, 24, 25, 26, 27]
        #    8  [23, 23, 23, 26, 27, 28, 29, 30]
        #                   | ----------|
        #    6  [24, 24, 24, 27, 28, 29, 30, 31]
        #    4  [25, 25, 25, 28, 29, 30, 31, 32]
        #                   | ----------|
        #    2  [...]

        #       Y=12 13  14                  19
        # X=14  [30, 30, 30, 33, 34, 35, 36, 37]
        #   12  [30, 30, 30, 33, 34, 35, 36, 37]
        #   10  [30, 30, 30, 33, 34, 35, 36, 37]
        #    8  [33, 33, 33, 36, 37, 38, 39, 40]
        #                   | ----------|
        #    6  [34, 34, 34, 37, 38, 39, 40, 41]
        #    4  [35, 35, 35, 38, 39, 40, 41, 42]
        #                   | ----------|
        # X= 2  [36, 36, 36, 39, 40, 41, 42, 43]

        #                           17       20
        expected_counts = [0]*17 + [1, 2, 2, 1, 0, 0, 0, 0, 0, 0,
                                    1, 2, 2, 1, 0, 0, 0, 0, 0, 0, 1, 2, 2, 1]
        #                          27 28 29 30                   37
        dvh = get_dvh(self.ss, self.dose, 1)
        diffl = dvh.differential
        # Counts are normalized to total, and to volume,
        # So undo that here for test dose grid.
        # 18=num dose voxels inside struct; 0.36=volume
        got_counts = diffl.counts * 18 / 0.36
        assert numpy.all(numpy.isclose(got_counts, expected_counts))

    def test_FF_decubitus_right(self):
        """Test DVH for feet-first decubitus right orientation."""
        self.dose.ImageOrientationPatient = [0, -1, 0, -1, 0, 0]
        self.dose.PixelSpacing = [2.0, 1.0]  # between Rows, Columns
        # original ipp = [2, 12, -20]
        # Then X = -r * dr + ipp[0], X decreases down the rows
        #  and Y = -c * dc + ipp[1], Y decreases across cols
        # (https://nipy.org/nibabel/dicom/dicom_orientation.html
        # #dicom-affine-formula)
        self.dose.ImagePositionPatient = [14, 19, 20]  # X Y Z top left
        # Below show contours box of (3, 14.5) - (7, 17.5) on dose grid
        #       Y=19 18  17  16  15  14  13  12
        # X=14  [10, 10, 10, 13, 14, 15, 16, 17],
        #   12  [10, 10, 10, 13, 14, 15, 16, 17]
        #   10  [10, 10, 10, 13, 14, 15, 16, 17]
        #    8  [13, 13, 13, 16, 17, 18, 19, 20]
        #               | ----------|
        #    6  [14, 14, 14, 17, 18, 19, 20, 21]
        #    4  [15, 15, 15, 18, 19, 20, 21, 22]
        #               | ----------|
        #    2  [16, 16, 16, 19, 20, 21, 22, 23]]

        #       Y=19 18  17  16  15  14  13  12
        # X=14  [20, 20, 20, 23, 24, 25, 26, 27]
        #   12  [20, 20, 20, 23, 24, 25, 26, 27]
        #   10  [20, 20, 20, 23, 24, 25, 26, 27]
        #    8  [23, 23, 23, 26, 27, 28, 29, 30]
        #               | ----------|
        #    6  [24, 24, 24, 27, 28, 29, 30, 31]
        #    4  [25, 25, 25, 28, 29, 30, 31, 32]
        #               | ----------|
        #    2  [...]

        #       Y=19 18  17  16  15  14  13  12
        # X=14  [30, 30, 30, 33, 34, 35, 36, 37]
        #   12  [30, 30, 30, 33, 34, 35, 36, 37]
        #   10  [30, 30, 30, 33, 34, 35, 36, 37]
        #    8  [33, 33, 33, 36, 37, 38, 39, 40]
        #               | ----------|
        #    6  [34, 34, 34, 37, 38, 39, 40, 41]
        #    4  [35, 35, 35, 38, 39, 40, 41, 42]
        #               | ----------|
        # X= 2  [36, 36, 36, 39, 40, 41, 42, 43]

        #                          14 15 16       19             24
        expected_counts = [0]*14 + [1, 1, 0, 1, 2, 1, 0, 0, 0, 0, 1, 1, 0, 1,
                                    2, 1, 0, 0, 0, 0, 1, 1, 0, 1, 2, 1]
        #                                            34
        dvh = get_dvh(self.ss, self.dose, 1)
        diffl = dvh.differential
        # Counts are normalized to total, and to volume,
        # So undo that here for test dose grid.
        # 18=num dose voxels inside struct; 0.36=volume
        got_counts = diffl.counts * 18 / 0.36
        assert numpy.all(numpy.isclose(got_counts, expected_counts))

    def test_FF_decubitus_right_structure_extents(self):
        """Test DVH for FF decubitus Rt orientation using structure extents."""
        self.dose.ImageOrientationPatient = [0, -1, 0, -1, 0, 0]
        self.dose.PixelSpacing = [2.0, 1.0]  # between Rows, Columns
        self.dose.ImagePositionPatient = [14, 19, 20]  # X Y Z top left
        # see grid from test_FF_decubitus_right
        #                          14 15 16       19             24
        expected_counts = [0]*14 + [1, 1, 0, 1, 2, 1, 0, 0, 0, 0, 1, 1, 0,
                                    1, 2, 1, 0, 0, 0, 0, 1, 1, 0, 1, 2, 1]
        #                                               34
        dvh = get_dvh(self.ss, self.dose, 1, use_structure_extents=True)
        diffl = dvh.differential
        # Counts are normalized to total, and to volume,
        # So undo that here for test dose grid.
        # 18=num dose voxels inside struct; 0.36=volume
        got_counts = diffl.counts * 18 / 0.36
        assert numpy.all(numpy.isclose(got_counts, expected_counts))

    def test_FF_decubitus_left(self):
        """Test DVH for feet-first decubitus left orientation."""
        self.dose.ImageOrientationPatient = [0, 1, 0, 1, 0, 0]
        self.dose.PixelSpacing = [2.0, 1.0]  # between Rows, Columns
        # original ipp = [2, 12, -20]
        # Then X = r * dr + ipp[0], X increases down the rows
        #  and Y = c * dc + ipp[1], Y increases across cols
        # (https://nipy.org/nibabel/dicom/dicom_orientation.html
        # #dicom-affine-formula)

        # In this test, we also shift Z so three structure planes use the
        #    first three dose planes rather than the middle three,
        #    just to ensure asymmetry in z direction is checked.
        #    Note, planes should really be reversed in pixel array, but doesn't
        #    matter since contour is identical on each slice.
        self.dose.ImagePositionPatient = [2, 12, 10]  # X Y Z top left
        # Below show contours box of (3, 14.5) - (7, 17.5) on dose grid
        #      Y=12  13  14  15  16  17      19
        # X=2   [ 0,  0,  0,  3,  4,  5,  6,  7],
        #                   |-----------|
        #   4   [ 0,  0,  0,  3,  4,  5,  6,  7]
        #   6   [ 0,  0,  0,  3,  4,  5,  6,  7]
        #                   |-----------|
        #   8   [ 3,  3,  3,  6,  7,  8,  9, 10]
        #  10   [ 4,  4,  4,  7,  8,  9, 10, 11]
        #  12   [ 5,  5,  5,  8,  9, 10, 11, 12]
        #  14   [ 6,  6,  6,  9, 10, 11, 12, 13]]

        #      Y=12  13  14                  19
        # X=2   [10, 10, 10, 13, 14, 15, 16, 17],
        #                   |-----------|
        #   4   [10, 10, 10, 13, 14, 15, 16, 17]
        #   6   [10, 10, 10, 13, 14, 15, 16, 17]
        #                   |-----------|
        #   8   [13, 13, 13, 16, 17, 18, 19, 20]
        #  10   [14, 14, 14, 17, 18, 19, 20, 21]
        #  12   [15, 15, 15, 18, 19, 20, 21, 22]
        #  14   [16, 16, 16, 19, 20, 21, 22, 23]]

        #      Y=12  13  14                  19
        # X=2   [20, 20, 20, 23, 24, 25, 26, 27]
        #                   |-----------|
        #   4   [20, 20, 20, 23, 24, 25, 26, 27]
        #   6   [20, 20, 20, 23, 24, 25, 26, 27]
        #                   |-----------|
        #   8   [23, 23, 23, 26, 27, 28, 29, 30]
        #  10   [24, 24, 24, 27, 28, 29, 30, 31]
        #  12   [25, 25, 25, 28, 29, 30, 31, 32]
        #  14   [...]

        #                          3
        expected_counts = [0]*3 + [2, 2, 2, 0, 0, 0, 0, 0, 0, 0,
                                   2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2]
        #                         13                            23
        dvh = get_dvh(self.ss, self.dose, 1)
        diffl = dvh.differential
        # Counts are normalized to total, and to volume,
        # So undo that here for test dose grid.
        # 18=num dose voxels inside struct; 0.36=volume
        got_counts = diffl.counts * 18 / 0.36
        assert numpy.all(numpy.isclose(got_counts, expected_counts))

    def test_empty_dose_grid(self):
        # See #274, prior to fixes this raised IndexError from
        #  get_interpolated_dose() getting empty array from GetDoseGrid()
        # Use z value to force no dose grid at that value
        #  Otherwise make like decub example
        self.dose.ImagePositionPatient = [2, 19, -1020]  # X Y Z top left
        self.dose.PixelSpacing = [2.0, 1.0]  # between Rows, Columns

        # 1 = roi number
        dvh = get_dvh(self.ss, self.dose, 1, use_structure_extents=True)
        self.assertTrue('Empty DVH' in dvh.notes)

    def test_not_implemented_orientations(self):
        """Test unhandled orientations raise NotImplementedError."""
        self.dose.ImageOrientationPatient = [0.7071, 0.7071, 0, 1, 0, 0]
        with self.assertRaises(NotImplementedError):
            _ = get_dvh(self.ss, self.dose, 1)


if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
