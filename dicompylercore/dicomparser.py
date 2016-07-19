#!/usr/bin/env python
# -*- coding: utf-8 -*-
# dicomparser.py
"""Class that parses and returns formatted DICOM RT data."""
# Copyright (c) 2009-2016 Aditya Panchal
# Copyright (c) 2009-2010 Roy Keyes
# This file is part of dicompyler-core, released under a BSD license.
#    See the file license.txt included with this distribution, also
#    available at https://github.com/dicompyler/dicompyler-core/


import logging
import numpy as np
try:
    from pydicom.dicomio import read_file
    from pydicom.dataset import Dataset
except ImportError:
    from dicom import read_file
    from dicom.dataset import Dataset
import random
from numbers import Number
from six import PY2, iterkeys
from six.moves import range
from dicompylercore import dvh, util

pil_available = True
try:
    from PIL import Image
except:
    pil_available = False

logger = logging.getLogger('dicompylercore.dicomparser')

class DicomParser:
    """Parses DICOM / DICOM RT files."""

    def __init__(self, dataset):

        if isinstance(dataset, Dataset):
            self.ds = dataset
        elif isinstance(dataset, str):
            try:
                self.ds = \
                    read_file(dataset, defer_size=100, force=True)
            except:
                # Raise the error for the calling method to handle
                raise
            else:
                # Sometimes DICOM files may not have headers,
                # but they should always have a SOPClassUID
                # to declare what type of file it is.
                # If the file doesn't have a SOPClassUID,
                # then it probably isn't DICOM.
                if not "SOPClassUID" in self.ds:
                    raise AttributeError
        else:
            raise AttributeError

######################## SOP Class and Instance Methods #######################

    def GetSOPClassUID(self):
        """Determine the SOP Class UID of the current file."""

        if (self.ds.SOPClassUID == '1.2.840.10008.5.1.4.1.1.481.2'):
            return 'rtdose'
        elif (self.ds.SOPClassUID == '1.2.840.10008.5.1.4.1.1.481.3'):
            return 'rtss'
        elif (self.ds.SOPClassUID == '1.2.840.10008.5.1.4.1.1.481.5'):
            return 'rtplan'
        elif (self.ds.SOPClassUID == '1.2.840.10008.5.1.4.1.1.2'):
            return 'ct'
        else:
            return None

    def GetSOPInstanceUID(self):
        """Determine the SOP Class UID of the current file."""

        return self.ds.SOPInstanceUID

    def GetStudyInfo(self):
        """Return the study information of the current file."""

        study = {}
        if 'StudyDescription' in self.ds:
            desc = self.ds.StudyDescription
        else:
            desc = 'No description'
        study['description'] = desc
        if 'StudyDate' in self.ds:
            date = self.ds.StudyDate
        else:
            date = None
        study['date'] = date
        # Don't assume that every dataset includes a study UID
        if 'StudyInstanceUID' in self.ds:
            study['id'] = self.ds.StudyInstanceUID
        else:
            study['id'] = str(random.randint(0, 65535))

        return study

    def GetSeriesInfo(self):
        """Return the series information of the current file."""

        series = {}
        if 'SeriesDescription' in self.ds:
            desc = self.ds.SeriesDescription
        else:
            desc = 'No description'
        series['description'] = desc
        series['id'] = self.ds.SeriesInstanceUID
        # Don't assume that every dataset includes a study UID
        series['study'] = self.ds.SeriesInstanceUID
        if 'StudyInstanceUID' in self.ds:
            series['study'] = self.ds.StudyInstanceUID
        series['referenceframe'] = self.ds.FrameOfReferenceUID \
            if 'FrameOfReferenceUID' in self.ds \
            else str(random.randint(0, 65535))
        if 'Modality' in self.ds:
            series['modality'] = self.ds.Modality
        else:
            series['modality'] = 'OT'

        return series

    def GetReferencedSeries(self):
        """Return the SOP Class UID of the referenced series."""

        if "ReferencedFrameOfReferenceSequence" in self.ds:
            frame = self.ds.ReferencedFrameOfReferenceSequence
            if "RTReferencedStudySequence" in frame[0]:
                study = frame[0].RTReferencedStudySequence[0]
                if "RTReferencedSeriesSequence" in study:
                    if "SeriesInstanceUID" in \
                            study.RTReferencedSeriesSequence[0]:
                        series = study.RTReferencedSeriesSequence[0]
                        return series.SeriesInstanceUID
        else:
            return ''

    def GetFrameOfReferenceUID(self):
        """Determine the Frame of Reference UID of the current file."""

        if 'FrameOfReferenceUID' in self.ds:
            return self.ds.FrameOfReferenceUID
        elif 'ReferencedFrameOfReferences' in self.ds:
            return self.ds.ReferencedFrameOfReferences[0].FrameOfReferenceUID
        else:
            return ''

    def GetReferencedStructureSet(self):
        """Return the SOP Class UID of the referenced structure set."""

        if "ReferencedStructureSetSequence" in self.ds:
            return self.ds.ReferencedStructureSetSequence[0].ReferencedSOPInstanceUID
        else:
            return ''

    def GetReferencedRTPlan(self):
        """Return the SOP Class UID of the referenced RT plan."""

        if "ReferencedRTPlanSequence" in self.ds:
            return self.ds.ReferencedRTPlanSequence[0].ReferencedSOPInstanceUID
        else:
            return ''

    def GetDemographics(self):
        """Return the patient demographics from a DICOM file."""

        # Set up some sensible defaults for demographics
        patient = {'name': 'None',
                   'id': 'None',
                   'birth_date': None,
                   'gender': 'Other'}
        if 'PatientName' in self.ds:
            if PY2:
                self.ds.decode()
            name = self.ds.PatientName
            patient['name'] = name
            patient['given_name'] = name.given_name
            patient['middle_name'] = name.middle_name
            patient['family_name'] = name.family_name
        if 'PatientID' in self.ds:
            patient['id'] = self.ds.PatientID
        if 'PatientSex' in self.ds:
            if (self.ds.PatientSex == 'M'):
                patient['gender'] = 'M'
            elif (self.ds.PatientSex == 'F'):
                patient['gender'] = 'F'
            else:
                patient['gender'] = 'O'
        if 'PatientBirthDate' in self.ds:
            if len(self.ds.PatientBirthDate):
                patient['birth_date'] = str(self.ds.PatientBirthDate)

        return patient


############################### Image Methods #################################

    def GetImageData(self):
        """Return the image data from a DICOM file."""

        data = {}

        if 'ImagePositionPatient' in self.ds:
            data['position'] = self.ds.ImagePositionPatient
        if 'ImageOrientationPatient' in self.ds:
            data['orientation'] = self.ds.ImageOrientationPatient
        if 'PixelSpacing' in self.ds:
            data['pixelspacing'] = self.ds.PixelSpacing
        else:
            data['pixelspacing'] = [1, 1]
        data['rows'] = self.ds.Rows
        data['columns'] = self.ds.Columns
        data['samplesperpixel'] = self.ds.SamplesPerPixel
        data['photometricinterpretation'] = self.ds.PhotometricInterpretation
        data['littlendian'] = \
            self.ds.file_meta.TransferSyntaxUID.is_little_endian
        if 'PatientPosition' in self.ds:
            data['patientposition'] = self.ds.PatientPosition
        data['frames'] = self.GetNumberOfFrames()

        return data

    def GetImageLocation(self):
        """Calculate the location of the current image slice."""

        ipp = self.ds.ImagePositionPatient
        iop = self.ds.ImageOrientationPatient

        normal = []
        normal.append(iop[1] * iop[5] - iop[2] * iop[4])
        normal.append(iop[2] * iop[3] - iop[0] * iop[5])
        normal.append(iop[0] * iop[4] - iop[1] * iop[3])

        loc = 0
        for i in range(0, len(normal)):
            loc += normal[i] * ipp[i]

        # The image location is inverted for Feet First images
        if 'PatientPosition' in self.ds:
            if ('ff' in self.ds.PatientPosition.lower()):
                loc = loc * -1

        return loc

    def GetImageOrientationType(self):
        """Get the orientation of the current image slice."""

        if 'ImageOrientationPatient' in self.ds:
            iop = np.array(self.ds.ImageOrientationPatient)

            orientations = [
                ["SA", np.array([1, 0, 0, 0, 1, 0])],      # supine axial
                ["PA", np.array([-1, 0, 0, 0, -1, 0])],    # prone axial
                ["SS", np.array([0, 1, 0, 0, 0, -1])],     # supine sagittal
                ["PS", np.array([0, -1, 0, 0, 0, -1])],    # prone sagittal
                ["SC", np.array([1, 0, 0, 0, 0, -1])],     # supine coronal
                ["PC", np.array([-1, 0, 0, 0, 0, -1])]     # prone coronal
            ]

            for o in orientations:
                if (not np.any(np.array(np.round(iop - o[1]), dtype=np.int32))):
                    return o[0]
        # Return N/A if the orientation was not found or could not be determined
        return "NA"

    def GetNumberOfFrames(self):
        """Return the number of frames in a DICOM image file."""

        frames = 1
        if 'NumberOfFrames' in self.ds:
            frames = self.ds.NumberOfFrames.real
        else:
            try:
                self.ds.pixel_array
            except:
                return 0
            else:
                if (self.ds.pixel_array.ndim > 2):
                    if (self.ds.SamplesPerPixel == 1) and not \
                       (self.ds.PhotometricInterpretation == 'RGB'):
                        frames = self.ds.pixel_array.shape[0]
        return frames

    def GetRescaleInterceptSlope(self):
        """Return the rescale intercept and slope if present."""

        intercept, slope = 0, 1
        if ('RescaleIntercept' in self.ds and 'RescaleSlope' in self.ds):
            intercept = self.ds.RescaleIntercept if \
                isinstance(self.ds.RescaleIntercept, Number) else 0
            slope = self.ds.RescaleSlope if \
                isinstance(self.ds.RescaleSlope, Number) else 1

        return intercept, slope

    def GetImage(self, window=0, level=0, size=None, background=False,
                 frames=0):
        """Return the image from a DICOM image storage file."""

        if not pil_available:
            print("Python imaging library not available." + \
                " Cannot generate images.")
            return

        # Return None if the Numpy pixel array cannot be accessed
        try:
            self.ds.pixel_array
        except:
            return Image.new('RGB', size, (0, 0, 0))

        # Samples per pixel are > 1 & RGB format
        if (self.ds.SamplesPerPixel > 1) and \
           (self.ds.PhotometricInterpretation == 'RGB'):

            # Little Endian
            if self.ds.file_meta.TransferSyntaxUID.is_little_endian:
                im = Image.frombuffer('RGB', (self.ds.Columns, self.ds.Rows),
                                      self.ds.PixelData, 'raw', 'RGB', 0, 1)
            # Big Endian
            else:
                im = Image.fromarray(np.rollaxis(
                    self.ds.pixel_array.transpose(), 0, 2))

        # Otherwise the image is monochrome
        else:
            if ((window == 0) and (level == 0)):
                window, level = self.GetDefaultImageWindowLevel()
            # Rescale the slope and intercept of the image if present
            intercept, slope = self.GetRescaleInterceptSlope()
            # Get the requested frame if multi-frame
            if (frames > 0):
                pixel_array = self.ds.pixel_array[frames]
            else:
                pixel_array = self.ds.pixel_array

            rescaled_image = pixel_array * slope + intercept

            image = self.GetLUTValue(rescaled_image, window, level)
            im = Image.fromarray(image).convert('L')

        # Resize the image if a size is provided
        if size:
            im.thumbnail(size, Image.ANTIALIAS)

        # Add a black background if requested
        if background:
            bg = Image.new('RGBA', size, (0, 0, 0, 255))
            bg.paste(im, ((size[0] - im.size[0]) / 2,
                     (size[1] - im.size[1]) / 2))
            return bg

        return im

    def GetDefaultImageWindowLevel(self):
        """Determine the default window/level for the DICOM image."""

        window, level = 0, 0
        if ('WindowWidth' in self.ds) and ('WindowCenter' in self.ds):
            if isinstance(self.ds.WindowWidth, float):
                window = self.ds.WindowWidth
            elif isinstance(self.ds.WindowWidth, list):
                if (len(self.ds.WindowWidth) > 1):
                    window = self.ds.WindowWidth[1]
            if isinstance(self.ds.WindowCenter, float):
                level = self.ds.WindowCenter
            elif isinstance(self.ds.WindowCenter, list):
                if (len(self.ds.WindowCenter) > 1):
                    level = self.ds.WindowCenter[1]

        if ((window, level) == (0, 0)):
            wmax = 0
            wmin = 0
            # Rescale the slope and intercept of the image if present
            intercept, slope = self.GetRescaleInterceptSlope()
            pixel_array = self.ds.pixel_array * slope + intercept

            if (pixel_array.max() > wmax):
                wmax = pixel_array.max()
            if (pixel_array.min() < wmin):
                wmin = pixel_array.min()
            # Default window is the range of the data array
            window = int(abs(wmax) + abs(wmin))
            # Default level is the range midpoint minus the window minimum
            level = int(window / 2 - abs(wmin))
        return window, level

    def GetLUTValue(self, data, window, level):
        """Apply the RGB Look-Up Table for the data and window/level value."""

        lutvalue = util.piecewise(data,
                                [data <= (level - 0.5 - (window - 1) / 2),
                                 data > (level - 0.5 + (window - 1) / 2)],
                                [0, 255, lambda data:
                                 ((data - (level - 0.5)) / (window-1) + 0.5) *
                                 (255 - 0)])
        # Convert the resultant array to an unsigned 8-bit array to create
        # an 8-bit grayscale LUT since the range is only from 0 to 255
        return np.array(lutvalue, dtype=np.uint8)

    def GetPatientToPixelLUT(self):
        """Get the image transformation matrix from the DICOM standard Part 3
            Section C.7.6.2.1.1"""

        di = self.ds.PixelSpacing[0]
        dj = self.ds.PixelSpacing[1]
        orientation = self.ds.ImageOrientationPatient
        position = self.ds.ImagePositionPatient

        m = np.matrix(
            [[orientation[0]*di, orientation[3]*dj, 0, position[0]],
            [orientation[1]*di, orientation[4]*dj, 0, position[1]],
            [orientation[2]*di, orientation[5]*dj, 0, position[2]],
            [0, 0, 0, 1]])

        x = []
        y = []
        for i in range(0, self.ds.Columns):
            imat = m * np.matrix([[i], [0], [0], [1]])
            x.append(float(imat[0]))
        for j in range(0, self.ds.Rows):
            jmat = m * np.matrix([[0], [j], [0], [1]])
            y.append(float(jmat[1]))

        return (np.array(x), np.array(y))

########################## RT Structure Set Methods ###########################

    def GetStructureInfo(self):
        """Return the patient demographics from a DICOM file."""

        structure = {}
        structure['label'] = self.ds.StructureSetLabel
        structure['date'] = self.ds.StructureSetDate
        structure['time'] = self.ds.StructureSetTime
        structure['numcontours'] = len(self.ds.ROIContourSequence)

        return structure

    def GetStructures(self):
        """Returns a dictionary of structures (ROIs)."""

        structures = {}

        # Determine whether this is RT Structure Set file
        if not (self.GetSOPClassUID() == 'rtss'):
            return structures

        # Locate the name and number of each ROI
        if 'StructureSetROISequence' in self.ds:
            for item in self.ds.StructureSetROISequence:
                data = {}
                number = int(item.ROINumber)
                data['id'] = number
                data['name'] = item.ROIName
                logger.debug("Found ROI #%s: %s", str(number), data['name'])
                structures[number] = data

        # Determine the type of each structure (PTV, organ, external, etc)
        if 'RTROIObservationsSequence' in self.ds:
            for item in self.ds.RTROIObservationsSequence:
                number = item.ReferencedROINumber
                structures[number]['type'] = item.RTROIInterpretedType

        # The coordinate data of each ROI is stored within ROIContourSequence
        if 'ROIContourSequence' in self.ds:
            for roi in self.ds.ROIContourSequence:
                number = roi.ReferencedROINumber

                # Generate a random color for the current ROI
                structures[number]['color'] = np.array((
                    random.randint(0, 255),
                    random.randint(0, 255),
                    random.randint(0, 255)), dtype=int)
                # Get the RGB color triplet for the current ROI if it exists
                if 'ROIDisplayColor' in roi:
                    # Make sure the color is not none
                    if not (roi.ROIDisplayColor is None):
                        color = roi.ROIDisplayColor
                    # Otherwise decode values separated by forward slashes
                    else:
                        value = roi[0x3006, 0x002a].repval
                        color = value.strip("'").split("/")
                    # Try to convert the detected value to a color triplet
                    try:
                        structures[number]['color'] = \
                            np.array(color, dtype=int)
                    # Otherwise fail and fallback on the random color
                    except:
                        logger.debug(
                            "Unable to decode display color for ROI #%s",
                            str(number))

                # Determine whether the ROI has any contours present
                if 'ContourSequence' in roi:
                    structures[number]['empty'] = False
                else:
                    structures[number]['empty'] = True

        return structures

    def GetStructureCoordinates(self, roi_number):
        """Get the list of coordinates for each plane of the structure."""

        planes = {}
        # The coordinate data of each ROI is stored within ROIContourSequence
        if 'ROIContourSequence' in self.ds:
            for roi in self.ds.ROIContourSequence:
                if (roi.ReferencedROINumber == int(roi_number)):
                    if 'ContourSequence' in roi:
                        # Locate the contour sequence for each referenced ROI
                        for c in roi.ContourSequence:
                            # For each plane, initialize a new plane dict
                            plane = {}

                            # Determine all the plane properties
                            plane['type'] = c.ContourGeometricType
                            plane['num_points'] = int(c.NumberOfContourPoints)
                            plane['data'] = \
                                self.GetContourPoints(c.ContourData)

                            # Each plane which coincides with an image slice
                            # will have a unique ID
                            if 'ContourImageSequence' in c:
                                # Add each plane to the planes dict
                                # of the current ROI
                                z = str(round(plane['data'][0][2], 2)) + '0'
                                if z not in planes:
                                    planes[z] = []
                                planes[z].append(plane)

        return planes

    def GetContourPoints(self, array):
        """Parses an array of xyz points & returns a array of point dicts."""

        n = 3
        return [array[i:i+n] for i in range(0, len(array), n)]

    def CalculatePlaneThickness(self, planesDict):
        """Calculates the plane thickness for each structure."""

        planes = []

        # Iterate over each plane in the structure
        for z in iterkeys(planesDict):
            planes.append(float(z))
        planes.sort()

        # Determine the thickness
        thickness = 10000
        for n in range(0, len(planes)):
            if (n > 0):
                newThickness = planes[n] - planes[n-1]
                if (newThickness < thickness):
                    thickness = newThickness

        # If the thickness was not detected, set it to 0
        if (thickness == 10000):
            thickness = 0

        return thickness

############################## RT Dose Methods ###############################

    def HasDVHs(self):
        """Returns whether dose-volume histograms (DVHs) exist."""

        if not "DVHSequence" in self.ds:
            return False
        else:
            return True

    def GetDVHs(self):
        """Returns cumulative dose-volume histograms (DVHs)."""

        self.dvhs = {}

        if self.HasDVHs():
            for item in self.ds.DVHSequence:
                # Make sure that the DVH has a referenced structure / ROI
                if not 'DVHReferencedROISequence' in item:
                    continue
                number = item.DVHReferencedROISequence[0].ReferencedROINumber
                logger.debug("Found DVH for ROI #%s", str(number))
                self.dvhs[number] = dvh.DVH.from_dicom_dvh(self.ds, number)

        return self.dvhs

    def GetDoseGrid(self, z=0, threshold=0.5):
        """
        Return the 2d dose grid for the given slice position (mm).

        :param z:           Slice position in mm.
        :param threshold:   Threshold in mm to determine the max difference
                            from z to the closest dose slice without
                            using interpolation.
        :return:            An numpy 2d array of dose points.
        """

        # If this is a multi-frame dose pixel array,
        # determine the offset for each frame
        if 'GridFrameOffsetVector' in self.ds:
            z = float(z)
            # Get the initial dose grid position (z) in patient coordinates
            imagepatpos = self.ds.ImagePositionPatient[2]
            orientation = self.ds.ImageOrientationPatient[0]
            # Add the position to the offset vector to determine the
            # z coordinate of each dose plane
            planes = orientation * np.array(self.ds.GridFrameOffsetVector) + \
                imagepatpos
            frame = -1
            # Check to see if the requested plane exists in the array
            if (np.amin(np.fabs(planes - z)) < threshold):
                frame = np.argmin(np.fabs(planes - z))
            # Return the requested dose plane, since it was found
            if not (frame == -1):
                return self.ds.pixel_array[frame]
            # Check if the requested plane is within the dose grid boundaries
            elif ((z < np.amin(planes)) or (z > np.amax(planes))):
                return []
            # The requested plane was not found, so interpolate between planes
            else:
                # Determine the upper and lower bounds
                umin = np.fabs(planes - z)
                ub = np.argmin(umin)
                lmin = umin.copy()
                # Change the min value to the max so we can find the 2nd min
                lmin[ub] = np.amax(umin)
                lb = np.argmin(lmin)
                # Fractional distance of dose plane between upper & lower bound
                fz = (z - planes[lb]) / (planes[ub] - planes[lb])
                plane = self.InterpolateDosePlanes(
                    self.ds.pixel_array[ub], self.ds.pixel_array[lb], fz)
                return plane
        else:
            return []

    def InterpolateDosePlanes(self, uplane, lplane, fz):
        """Interpolates a dose plane between two bounding planes at the given
           relative location.

        :param uplane:      Upper dose plane boundary.
        :param uplane:      Lower dose plane boundary.
        :param fz:          Fractional distance from the bottom to the top,
                            where the new plane is located.
                            E.g. if fz = 1, the plane is at the upper plane,
                            fz = 0, it is at the lower plane.
        :return:            An numpy 2d array of the interpolated dose plane.
        """

        # A simple linear interpolation
        doseplane = fz*uplane + (1.0 - fz)*lplane

        return doseplane

    def GetIsodosePoints(self, z=0, level=100, threshold=0.5):
        """Return points for the given isodose level and slice position
           from the dose grid.

        :param z:           Slice position in mm.
        :param threshold:   Threshold in mm to determine the max difference
                            from z to the closest dose slice without
                            using interpolation.
        :param level:       Isodose level in scaled form
                            (multiplied by self.ds.DoseGridScaling)
        :return:            An array of tuples representing isodose points.
        """

        plane = self.GetDoseGrid(z, threshold)
        isodose = (plane >= level).nonzero()
        return list(zip(isodose[1].tolist(), isodose[0].tolist()))

    def GetDoseData(self):
        """Return the dose data from a DICOM RT Dose file."""

        data = self.GetImageData()
        data['doseunits'] = self.ds.DoseUnits
        data['dosetype'] = self.ds.DoseType
        data['dosecomment'] = self.ds.DoseComment \
            if 'DoseComment' in self.ds else ''
        data['dosesummationtype'] = self.ds.DoseSummationType
        data['dosegridscaling'] = self.ds.DoseGridScaling
        data['dosemax'] = float(self.ds.pixel_array.max())
        data['lut'] = self.GetPatientToPixelLUT()
        data['fraction'] = ''
        if "ReferencedRTPlanSequence" in self.ds:
            plan = self.ds.ReferencedRTPlanSequence[0]
            if "ReferencedFractionGroupSequence" in plan:
                data['fraction'] = \
                plan.ReferencedFractionGroupSequence[0].ReferencedFractionGroupNumber

        return data

    def GetReferencedBeamNumber(self):
        """Return the referenced beam number (if it exists) from RT Dose."""

        beam = None
        if "ReferencedRTPlanSequence" in self.ds:
            rp = self.ds.ReferencedRTPlanSequence[0]
            if "ReferencedFractionGroupSequence" in rp:
                rf = rp.ReferencedFractionGroupSequence[0]
                if "ReferencedBeamSequence" in rf:
                    if "ReferencedBeamNumber" in rf.ReferencedBeamSequence[0]:
                        beam = \
                            rf.ReferencedBeamSequence[0].ReferencedBeamNumber

        return beam

############################## RT Plan Methods ###############################

    def GetPlan(self):
        """Returns the plan information."""

        self.plan = {}

        self.plan['label'] = self.ds.RTPlanLabel
        self.plan['date'] = self.ds.RTPlanDate
        self.plan['time'] = self.ds.RTPlanTime
        self.plan['name'] = ''
        self.plan['rxdose'] = 0
        if "DoseReferenceSequence" in self.ds:
            for item in self.ds.DoseReferenceSequence:
                if item.DoseReferenceStructureType == 'SITE':
                    self.plan['name'] = "N/A"
                    if "DoseReferenceDescription" in item:
                        self.plan['name'] = item.DoseReferenceDescription
                    if 'TargetPrescriptionDose' in item:
                        rxdose = item.TargetPrescriptionDose * 100
                        if (rxdose > self.plan['rxdose']):
                            self.plan['rxdose'] = rxdose
                elif item.DoseReferenceStructureType == 'VOLUME':
                    if 'TargetPrescriptionDose' in item:
                        self.plan['rxdose'] = item.TargetPrescriptionDose * 100
        if (("FractionGroupSequence" in self.ds) and (self.plan['rxdose'] == 0)):
            fg = self.ds.FractionGroupSequence[0]
            if ("ReferencedBeams" in fg) and \
               ("NumberofFractionsPlanned" in fg):
                beams = fg.ReferencedBeams
                fx = fg.NumberofFractionsPlanned
                for beam in beams:
                    if "BeamDose" in beam:
                        self.plan['rxdose'] += beam.BeamDose * fx * 100
        self.plan['rxdose'] = int(self.plan['rxdose'])
        self.plan['brachy'] = False
        if ("BrachyTreatmentTechnique" in self.ds) or \
                ("BrachyTreatmentType" in self.ds):
            self.plan['brachy'] = True
        return self.plan

    def GetReferencedBeamsInFraction(self, fx=0):
        """Return the referenced beams from the specified fraction."""

        beams = {}
        if ("BeamSequence" in self.ds):
            bdict = self.ds.BeamSequence
        elif ("IonBeamSequence" in self.ds):
            bdict = self.ds.IonBeamSequence
        else:
            return beams

        # Obtain the beam information
        for b in bdict:
            beam = {}
            beam['name'] = b.BeamName if "BeamName" in b else ""
            beam['description'] = b.BeamDescription \
                if "BeamDescription" in b else ""
            beams[b.BeamNumber.real] = beam

        # Obtain the referenced beam info from the fraction info
        if ("FractionGroupSequence" in self.ds):
            fg = self.ds.FractionGroupSequence[fx]
            if ("ReferencedBeamSequence" in fg):
                rb = fg.ReferencedBeamSequence
                nfx = fg.NumberOfFractionsPlanned
                for b in rb:
                    if "BeamDose" in b:
                        beams[b.ReferencedBeamNumber]['dose'] = \
                            b.BeamDose * nfx * 100
        return beams
