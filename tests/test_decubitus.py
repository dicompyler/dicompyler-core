import unittest
import pydicom
from pydicom.dataset import Dataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import generate_uid, PYDICOM_ROOT_UID
import pydicom.uid
import numpy
from dicompylercore.dvhcalc import get_dvh

"""
Set up dose grid like:
    2 mm x 2 mm pixels, 10x10 grid
    In patient coords
    across (cols, X) from 1 to 21, center of pixels from 2 to 20
    down (rows, Y) from 11 to 31, centers from 12 to 30
    slices -20, -10, 0, 10, 20

    So ... ImagePositionPatient = (2, 12, -20),
    then slice offsets (0, 10, 20, 30, 40)

    Then make contour on the middle three slices, with rectangle coords:
    (3, 17), (7, 17), (7, 23), (3, 23)    (2 voxels across x, 3 down y)


"""

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
CT_iop = [0, -1, 0, 1, 0, 0]  # Feet-first decub left
#   [ 0, -1,  0, -1,  0,  0] # FF dec rt    # [1, 0, 0,  0, 1, 0]  Standard
DOSE_ipp = CT_ipp
DOSE_iop = CT_iop

def basic_file_meta(class_UID):
    # File meta info data elements
    file_meta = FileMetaDataset()
    file_meta.FileMetaInformationVersion = b'\x00\x01'
    file_meta.MediaStorageSOPClassUID = class_UID
    file_meta.MediaStorageSOPInstanceUID = generate_uid()
    file_meta.TransferSyntaxUID = '1.2.840.10008.1.2'
    file_meta.ImplementationClassUID = PYDICOM_ROOT_UID
    return file_meta


def fake_rtdose():
    # Main data elements
    file_meta = basic_file_meta(pydicom.uid.RTDoseStorage)
    ds = Dataset()
    ds.SOPClassUID = pydicom.uid.RTDoseStorage
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
    ds.Rows = 10
    ds.Columns = 10
    ds.PixelSpacing = [2.0, 2.0]
    ds.BitsAllocated = 32
    ds.BitsStored = 32
    ds.HighBit = 31
    ds.PixelRepresentation = 0
    ds.DoseUnits = 'GY'
    ds.DoseType = 'PHYSICAL'
    ds.DoseSummationType = 'PLAN'
    ds.GridFrameOffsetVector = [0.0, 10.0, 20.0, 30.0, 40.0]
    ds.DoseGridScaling = '1.0'

    # Referenced RT Plan Sequence
    refd_rt_plan_sequence = Sequence()
    ds.ReferencedRTPlanSequence = refd_rt_plan_sequence

    # Referenced RT Plan Sequence: Referenced RT Plan 1
    refd_rt_plan1 = Dataset()
    refd_rt_plan1.ReferencedSOPClassUID = '1.2.840.10008.5.1.4.1.1.481.5'
    refd_rt_plan1.ReferencedSOPInstanceUID = PLAN_iUID
    refd_rt_plan_sequence.append(refd_rt_plan1)


    arr = numpy.zeros((ds.NumberOfFrames, ds.Rows, ds.Columns), dtype=numpy.uint32)

    for slice in range(len(arr)):
        for row in range(len(arr[0])):
            for col in range(len(arr[0][0])):
                arr[slice, row, col] = slice*10 + (row > 5)*row + (col > 5) * (col)
    ds.PixelData = arr.tobytes()

    ds.file_meta = file_meta
    ds.is_implicit_VR = True
    ds.is_little_endian = True

    return ds


def fake_ss():
    file_meta = basic_file_meta(pydicom.uid.RTStructureSetStorage)
    # Main data elements
    ds = Dataset()
    ds.SpecificCharacterSet = 'ISO_IR 192'
    ds.InstanceCreationDate = '20220101'
    ds.SOPClassUID = pydicom.uid.RTStructureSetStorage
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
    contour_image1.ReferencedSOPClassUID = pydicom.uid.CTImageStorage
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

    for i in range(1, NUM_SLICES-1):  # contours on all but first and last slice
        contour = Dataset()

        # # Contour Image Sequence
        # contour_image_sequence = Sequence()
        # contour.ContourImageSequence = contour_image_sequence

        # # Contour Image Sequence: Contour Image 1
        # contour_image1 = Dataset()
        # contour_image1.ReferencedSOPClassUID = '1.2.840.10008.5.1.4.1.1.2'
        # contour_image1.ReferencedSOPInstanceUID = '1.3.6.1.4.1.22361.140839833039574.112288227.1641236932544.0'
        # contour_image_sequence.append(contour_image1)

        contour.ContourGeometricType = 'CLOSED_PLANAR'

        z = SLICE_Z[i]
        contour.ContourData = [3, 17, z,  7, 17, z,  7, 23, z,  3, 23, z]
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

    def test_decub(self):
        """Test that DVH is calculated correctly for decubitus orientation"""
        dvh = get_dvh(self.ss, self.dose, 1)
        self.assertTrue('Empty DVH' not in dvh.notes)
        # XXX add specific DVH results in

    def test_empty_dose_grid(self):
        # See #274, prior to fixes this raised IndexError from
        #  get_interpolated_dose() getting empty array from GetDoseGrid()
        self.dose.ImagePositionPatient[2] -= 1000 # force no dose at struct z
        dvh = get_dvh(self.ss, self.dose, 1, use_structure_extents=True)  # 1 = roi number
        self.assertTrue('Empty DVH' in dvh.notes)

    def test_not_implemented_orientations(self):
        """Test unhandled orientations raise NotImplementedError"""
        self.dose.ImageOrientationPatient = [0.7071, 0.7071, 0, 1, 0, 0]
        with self.assertRaises(NotImplementedError):
            dvh = get_dvh(self.ss, self.dose, 1)



if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
