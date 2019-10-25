#!/usr/bin/env python3

"""Plots data associated with Elastica simulations

Used to reproduce various figures and renderings (in-part)
from the paper:
Modeling and simulation of complex dynamic musculoskeletal architectures
Nature Communications, 2019
"""
__license__ = "MIT License, see LICENSE for details"
__copyright__ = "Copyright (C) 2019 MattiaLab"

# System imports
import argparse
import os
import sys
from itertools import tee

import matplotlib.style as mplstyle  #
import numpy as np  #
from matplotlib import pyplot as plt  #
from matplotlib.colors import to_rgb
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.mplot3d.art3d import Line3DCollection  #
from scipy.linalg import norm  #

# set the backend first
# import matplotlib # isort:skip
# matplotlib.use("TkAgg") # isort:skip
# from matplotlib.collections import LineCollection #


# Turned off because slow
# plt.rcParams['text.usetex'] = 'True'
# plt.rcParams['font.serif'] = 'Fira Sans'
plt.rcParams["font.size"] = 14
plt.rcParams["axes.labelsize"] = 14
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titlesize"] = 16
plt.rcParams["xtick.labelsize"] = 12
plt.rcParams["ytick.labelsize"] = 12
mplstyle.use("seaborn-whitegrid")
# plt.rc('grid', color='#397939', linewidth=1, linestyle='--')


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    first, second = tee(iterable)
    next(second, None)
    return zip(first, second)


def set_axes_equal(axis, data=None):
    """Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
    ax: a matplotlib axis, e.g., as output from plt.gca().

    Credits: https://stackoverflow.com/a/50664367
    """
    if data is None:
        limits = np.array([axis.get_xlim3d(), axis.get_ylim3d(), axis.get_zlim3d()])
    else:
        limits = np.array(
            [
                [data[0, :].min(), data[0, :].max()],
                [data[1, :].min(), data[1, :].max()],
                [data[2, :].min(), data[2, :].max()],
            ]
        )

    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    axis.set_xlim3d([origin[0] - radius, origin[0] + radius])
    axis.set_ylim3d([origin[1] - radius, origin[1] + radius])
    axis.set_zlim3d([origin[2] - radius, origin[2] + radius])


def detect_files(input_folder, prefix, suffix):
    """ Detect all files in the input folder
    """
    if not os.path.isdir(input_folder):
        raise FileNotFoundError("{} is not a valid folder".format(input_folder))

    # Finds all files input_folder/prefix*.suffix
    import re

    # matches prefix__ (bignumber).suffix
    exp = r"{prefix}[\s_]*(\d*){meta}{suffix}".format(
        prefix=prefix, meta="\\", suffix=suffix
    )
    compiled_regexp = re.compile(exp)

    matched_files = [f for f in os.listdir(input_folder) if compiled_regexp.search(f)]
    matched_files.sort()

    # Extract the filenumber as an integer
    try:
        file_tags = [int(compiled_regexp.match(f).group(1)) for f in matched_files]
    except ValueError:
        file_tags = []

    # Put the full address
    matched_files = [os.path.join(input_folder, f) for f in matched_files]

    # returns matched files and matched file numbers as strings
    return matched_files, file_tags


class Arrow3D(FancyArrowPatch):
    """3D Arrow plot for drawing vector, from https://stackoverflow.com/a/22867877"""

    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)


class CustomFormatter(
    argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter
):
    pass


class FigProperties:  # pylint: disable=R0903
    width = 900
    height = 600
    dpi = 100

    @staticmethod
    def figsize():
        return (
            FigProperties.width / float(FigProperties.dpi),
            FigProperties.height / float(FigProperties.dpi),
        )


class DummyPlotter:
    """ Placeholder when we need to skip plotting
    """

    # pylint: disable=R0913
    def __init__(
        self, input_folder, output_folder, save_file, force_flag, display_flag
    ):
        pass

    def process(self):
        pass

    def plot(self, axis, data, color=(31 / 255, 119 / 255, 180 / 255)):
        pass

    def animate(self):
        pass


class ThreeDimensionalPlotter:
    """ Base class for all three-dimensional plotting
    """

    def __init__(
        self, input_folder, output_folder, save_file, force_flag, display_flag
    ):

        """ Figure attributes """
        self.fig = plt.figure(figsize=FigProperties.figsize())

        self.ax = self.fig.add_subplot(111, projection="3d")
        self.ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        self.ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        self.ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        self.ax.grid(True)
        self.ax.set_xlabel("X")
        self.ax.set_ylabel("Y")
        self.ax.set_zlabel("Z")

        """ Save file attributes """
        self.output_folder = output_folder
        self.savefile_name, self.savefile_ext = os.path.splitext(save_file)
        self.display_flag = display_flag

        # Assume its a png if no extension given by default
        if not self.savefile_ext:
            self.savefile_ext = ".png"

        """ Source file attributes and search """
        # Guaranteed to be sorted
        ntypes_files = len(self.file_metadata)
        # Get a collection (list of lists) to store all the source files
        src_file_collection = [[] for i in range(ntypes_files)]
        src_filetag_collection = [[] for i in range(ntypes_files)]

        for index, (prefix, suffix, _) in enumerate(self.file_metadata):
            src_file_collection[index], src_filetag_collection[index] = detect_files(
                input_folder, prefix, suffix
            )

        for i, j in pairwise(range(ntypes_files)):
            # These lists are ordered and the default comparison should be fine
            assert (
                src_filetag_collection[i] == src_filetag_collection[j]
            ), "Numbers don't match up"

        # All collectinos are the same, pop only the last one
        src_filetags = src_filetag_collection.pop()

        tgt_files, tgt_filetags = detect_files(
            self.output_folder, self.savefile_name, self.savefile_ext
        )
        # print(input_folder, src_files, src_filetags)
        # print(output_folder, tgt_files, tgt_filetags)

        # If uses forces to clear all images, process all files
        if force_flag:
            self.files_to_be_processed = src_file_collection
            self.filetags_to_be_processed = src_filetags

        # Else, rewrite files that only need to be updated
        else:
            if len(tgt_filetags) > len(src_filetags):
                # if more targets already, then there's something fishy
                # act as if the force_flag is set
                self.files_to_be_processed = src_file_collection
                self.filetags_to_be_processed = src_filetags
            else:
                # Calculate difference between the filetags in src and tgt
                # if tgt is not empty

                # # Get indices of sort list, see https://stackoverflow.com/a/6423325
                # sort_src_indices = sorted(range(len(src_filetags))
                #                           ,key=src_filetags.__getitem__)
                # sort_srctags = [src_filetags[i] for i in sort_src_indices]
                # sort_tgttags = sorted(tgt_filetags)

                def index(a, x):
                    """Locate the leftmost value exactly equal to x in a sorted list
                    See https://docs.python.org/3.7/library/bisect.html
                    """
                    import bisect

                    i = bisect.bisect_left(a, x)
                    if i != len(a) and a[i] == x:
                        return i
                    return -1

                # If tgttags is empty, default to -1
                # else what is the next element beyond the last target tag?
                src_index = (
                    -1 if not tgt_filetags else index(src_filetags, tgt_filetags[-1])
                )
                # print(src_index)

                self.files_to_be_processed = [
                    sublist[src_index + 1 :] for sublist in src_file_collection
                ]
                self.filetags_to_be_processed = src_filetags[src_index + 1 :]

        # print(self.files_to_be_processed, self.filetags_to_be_processed)

    def process(self):
        """ Loads all data, plots them and stores them into appropriately
        named figures.
        """

        # Load color metadata first as list of rgb tuples
        colors = [color for (_, _, color) in self.file_metadata]

        if self.display_flag:
            # Turn on interactive mode to persist figure
            plt.ion()
            # Show figure after persist
            plt.show()

        # From list of list (f) and list (g), get [f[i][0],f[i][1]] and g[i]
        for *src_file_names, src_file_tag in zip(
            *(self.files_to_be_processed), self.filetags_to_be_processed
        ):
            for i_seq, (src_file, color) in enumerate(zip(src_file_names, colors)):
                # Defaults loads to (ndata, 4) rather than (4,ndata)
                data = np.loadtxt(src_file).T
                self.plot(self.ax, data, color, i_seq)

            # self.ax.set_title("test at {}".format(src_file_tag))
            self.fig.canvas.draw()

            if self.display_flag:
                plt.pause(0.001)

            # last column is source files
            # src_file_tag = src_file_info[-1]

            filename = "{name}_{tag:05d}{ext}".format(
                name=self.savefile_name, tag=src_file_tag, ext=self.savefile_ext
            )
            # tight bounding box here screws up the video
            self.fig.savefig(
                os.path.join(self.output_folder, filename), dpi=FigProperties.dpi
            )

        # # Show only the last drawn figure to the user
        # plt.show()

    def plot(self, ax, data, color=(31 / 255, 119 / 255, 180 / 255), i_seq=None):
        """ Plots 3D data
        """
        # https://stackoverflow.com/a/34486703
        # Plot centerline
        if len(ax.lines) > i_seq:
            ax.lines[i_seq]._verts3d = data[:3, :]  # 0,1,2
        else:
            ax.plot(
                data[0, :],
                data[1, :],
                data[2, :],
                color=color,
                marker="o",
                markersize=data[3, 0] * 0.2,  # set marker based on radius
                linewidth=2.0,
            )

        pts_data, wire_data = calculate_cylinders(data)

        # Now plot wireframes
        # Each one adds two, so sequencing should reflect that
        if len(ax.collections) > 2 * i_seq:
            # First do pts_data
            ax.collections[2 * i_seq]._segments3d = convert_to_collection(pts_data)
            # Then do wire_data
            ax.collections[2 * i_seq + 1]._segments3d = convert_to_collection(wire_data)
        else:
            # first instance
            #                           transparency
            metadata = {"colors": color + (0.5,), "linewidth": 1.0}
            # ax.add_collection3d(convert_to_collection(pts_data, **metadata))
            ax.add_collection3d(
                Line3DCollection(convert_to_collection(pts_data), **metadata)
            )

            #                           transparency
            metadata = {"colors": color + (0.3,), "linewidth": 0.5}
            ax.add_collection3d(
                Line3DCollection(convert_to_collection(wire_data), **metadata)
            )

    # def plot(self, ax, data, color=(31 / 255, 119 / 255, 180 / 255), i_seq=None):
    #     """ Plots 3D data
    #     """
    #     # https://stackoverflow.com/a/34486703
    #     # Plot centerline
    #     if ax.lines:
    #         for line in ax.lines:
    #             line._verts3d = data[:-1, :]
    #     else:
    #         ax.plot(
    #             data[0, :],
    #             data[1, :],
    #             data[2, :],
    #             color=color,
    #             marker="o",
    #             linewidth=3.0,
    #         )

    #     pts_data, wire_data = calculate_cylinders(data)

    #     # Now plot wireframes
    #     if ax.collections:
    #         # First do pts_data
    #         ax.collections[0]._segments3d = convert_to_collection(pts_data)
    #         # Then do wire_data
    #         ax.collections[1]._segments3d = convert_to_collection(wire_data)
    #     else:
    #         # first instance
    #         #                           transparency
    #         metadata = {"colors": color + (1.0,), "linewidth": 2.0}
    #         # ax.add_collection3d(convert_to_collection(pts_data, **metadata))
    #         ax.add_collection3d(
    #             Line3DCollection(convert_to_collection(pts_data), **metadata)
    #         )

    #         # # Add data to scatter
    #         # ax.scatter(data[:, 0],
    #         #            data[:, 1],
    #         #            data[:, 2],
    #         #            color=colors[0],
    #         #            marker="o")

    #         #                           transparency
    #         metadata = {"colors": color + (0.5,), "linewidth": 1.5}
    #         ax.add_collection3d(
    #             Line3DCollection(convert_to_collection(wire_data), **metadata)
    #         )

    def animate(self):
        """ Animates using ffmpeg the figures created by savefig
        """
        import subprocess

        # Test if ffmpeg present, else stop processing
        ret_code = subprocess.call(["which", "ffmpeg"])

        if ret_code != 0:
            raise OSError("ffmpeg not found. Aborting now.")
        else:
            filenames = "{frame}_%05d{ext}".format(
                frame=os.path.join(self.output_folder, self.savefile_name),
                ext=self.savefile_ext,
            )

            # ffmpeg has its own glob matching facility
            subprocess.call(
                [
                    "ffmpeg",
                    "-y",  # overwrites files
                    "-i",
                    filenames,
                    "-framerate",
                    "30",
                    "-crf",
                    "24",
                    "-pix_fmt",
                    "yuv420p",
                    "output_video.mp4",
                ]
            )


def convert_to_collection(in_data, **kwargs):
    """ Converts a (3,*) np array into a MPL LineCollection
    """
    ndims = in_data.shape[0]
    if ndims == 3:
        points = np.array([in_data[0], in_data[1], in_data[2]]).T.reshape(-1, 1, 3)
        segs = np.concatenate([points[:-1], points[1:]], axis=1)
        # return Line3DCollection(segs, **kwargs)
    elif ndims == 2:
        points = np.array([in_data[0], in_data[1]]).T.reshape(-1, 1, 2)
        segs = np.concatenate([points[:-1], points[1:]], axis=1)
        # return LineCollection(segs, **kwargs)
    else:
        raise IndexError("Dimensions incorrect!")

    return segs


def calculate_cylinders(in_data):
    """ Calculates the cylinder coordinates given the centerline

    in_data : (n_dim + 1, N) in size, last dimension for radius
    """

    # print(in_data.shape)
    # Split dimensions and radius
    data = in_data[:-2, :]
    # Last one is time now discounted
    radius_data = in_data[-2, :-1]
    n_dim = data.shape[0]

    # Governing parameters
    n_pts = data.shape[1]
    n_axial = 2
    n_theta = 20
    n_wireframe_skip = 5
    n_wireframe_skip = n_wireframe_skip if n_wireframe_skip < n_theta else 1

    """ 1. Calculate tangent """
    tan = np.diff(data)

    # normalize
    mag_tan = norm(tan, ord=2, axis=0)
    tan /= mag_tan

    """ 2. Calculate normal and binormal """
    # Guess for binormal, fair enough to assume in z
    binormal = np.array([0.0, 0.5, 0.5])
    binormal /= norm(binormal)

    # prepare to broadcast it to 2D
    binormal = binormal[:, np.newaxis]

    # make vector perpendicular to v
    normal = np.cross(tan, np.tile(binormal, (1, n_pts - 1)), axisa=0, axisb=0, axisc=0)
    # normalize
    mag_n = norm(normal, ord=2, axis=0)
    normal /= mag_n

    # make unit vector perpendicular to v and n1
    binormal = np.cross(tan, normal, axisa=0, axisb=0, axisc=0)

    # Stack normal and binormal together for (2,ndim,N-1)
    directors = np.vstack((binormal[np.newaxis, :, :], normal[np.newaxis, :, :]))

    """ 3. Parametrize angles and centerline """
    # surface ranges over t from 0 to length of axis
    caxis = np.linspace(0.0, 1.0, n_axial)
    #  polar angle varies from 0 to 2*pi
    theta = np.linspace(0, 2 * np.pi, n_theta)

    """ 4. Direct them"""
    # Idea here is to direct the centerline in the tangent direction
    # and direct the cross section in the norm-binorm direction
    # and then take a linear combination of both at every cross section

    # scale t up according to mag_tan, to give (n_axial, N-1) array
    # t = t.reshape(n_axial, 1) * mag_tan.reshape(1, -1)
    caxis = np.einsum("i,j->ij", caxis, mag_tan)
    # Multiply by tangent to give a (ndim, n_axial, N-1) array, two ways
    # t = tan[:, np.newaxis, :] * t[np.newaxis, :, :]
    caxis = np.einsum("ik,jk->ijk", tan, caxis)

    # calculate cossin first to give a (2, n_theta) array
    cs_theta = np.vstack((np.cos(theta), np.sin(theta)))
    # scale cs_theta up to give (2, n_theta, N-1) array
    cs_theta = np.tile(cs_theta[:, :, np.newaxis], (1, 1, n_pts - 1))
    # multiply by elemental cross section radius, retain shape
    cs_theta = np.einsum("ijk,k->ijk", cs_theta, radius_data)
    # Multiply rcos and rsin by the appropriate normal, binorm
    # at every cross section & give (ndim, n_theta, N-1) array
    inplane = np.einsum("ijk,ilk->jlk", directors, cs_theta)

    # Get coordinate data for all points on perimeter
    # according to x^p_i = x^node_i + t*tan_i + r*cos*bn_i + r*sin*n_i
    # data converted from 1 to n-1 as there are only (n-1) elements
    # Gives a (ndim, n_axial, n_theta, N-1) array
    perimeter_pts = (
        data[:, np.newaxis, np.newaxis, :-1]
        + caxis[:, :, np.newaxis, :]
        + inplane[:, np.newaxis, :, :]
    )
    # Gives only select angular lines running from start circle to end circle
    wireframe_pts = perimeter_pts[:, [0, -1], ::n_wireframe_skip, :]
    # Note : here it is important to reshape perimeter pts to have
    # a sense of continuity in the azimuthal position
    # this is not important for the wireframe however
    return perimeter_pts.reshape(n_dim, n_axial * n_theta, -1), wireframe_pts


class SphericalJointPlotter(ThreeDimensionalPlotter):
    """ Plots the trajectory of the slithering sphericalJoint
    as simulation progresses
    """

    def __init__(
        self, input_folder, output_folder, save_file, force_flag, display_flag
    ):

        # # Body and wireframe color
        # self.file_metadata = [
        #     ("rod1", ".txt", (31 / 255, 119 / 255, 180 / 255)),
        #     ("rod2", ".txt", (200 / 255, 0 / 255, 0 / 255)),
        # ]

        # Body and wireframe color
        self.file_metadata = [
            ("rod1", ".txt", to_rgb("xkcd:bluish")),
            ("rod2", ".txt", to_rgb("xkcd:reddish")),
        ]

        super(SphericalJointPlotter, self).__init__(
            input_folder, output_folder, save_file, force_flag, display_flag
        )

        # Make any other changes to the figure here
        # Data-dependent maybe
        self.ax.set_xlim(-200, 200)
        self.ax.set_ylim(-150, 250)
        self.ax.set_zlim(00, 400)
        # self.ax.set_ylim(-0.125, 0.125)
        # self.ax.set_zlim(-0.2, 0.2)
        set_axes_equal(self.ax)

    def plot(self, ax, data, color=(31 / 255, 119 / 255, 180 / 255), i_seq=None):
        super(SphericalJointPlotter, self).plot(ax, data, color, i_seq)
        last_pos = data[0:3, -1]  # position
        time = data[4, 0]  # time
        scaling = 60.0  # purely for visualization

        # Below routine for plotting force as a line, and not arrow

        # # After the last possible dynamic update plot
        # if i_seq == 1:  # ie. i_seq == len(self.metadata)-1
        #     if len(ax.lines) > i_seq + 1:  # 2 * i_seq + 2
        #         if time > 0.2:
        #             f = 1.0 * np.array(
        #                 [
        #                     np.cos(0.5 * np.pi * (time - 0.2)),
        #                     0.0,
        #                     np.sin(0.5 * np.pi * (time - 0.2)),
        #                 ]
        #             )
        #         else: # time < 0.2, only update in position is needed
        #              f = np.array([0.0, 0.0, -2])

        #         f_with_pos = np.vstack((last_pos, last_pos + scaling * f)).T
        #         ax.lines[i_seq + 1]._verts3d = f_with_pos
        #     else:
        #         f = scaling * np.array([0.0, 0.0, -2])
        #         f_with_pos = np.vstack((last_pos, last_pos + f)).T
        #         ax.plot(
        #             f_with_pos[0, :],
        #             f_with_pos[1, :],
        #             f_with_pos[2, :],
        #             color="k",
        #             linewidth=1,
        #         )

        # After the last possible dynamic update plot
        if i_seq == 1:  # ie. i_seq == len(self.metadata)-1
            if ax.artists:
                if time > 0.2:
                    f = 1.0 * np.array(
                        [
                            np.cos(0.5 * np.pi * (time - 0.2)),
                            0.0,
                            np.sin(0.5 * np.pi * (time - 0.2)),
                        ]
                    )
                else:  # time < 0.2, only update in position is needed
                    f = np.array([0.0, 0.0, -2])
                arrow_end_pos = last_pos + scaling * f

                ax.artists[0]._verts3d = np.vstack((last_pos, arrow_end_pos)).T
            else:
                f = np.array([0.0, 0.0, -2])
                arrow_end_pos = last_pos + scaling * f
                a = Arrow3D(
                    [last_pos[0], arrow_end_pos[0]],
                    [last_pos[1], arrow_end_pos[1]],
                    [last_pos[2], arrow_end_pos[2]],
                    mutation_scale=10,
                    lw=1,
                    arrowstyle="-|>",
                    color="k",
                )
                ax.add_artist(a)


class HingeJointPlotter(SphericalJointPlotter):
    """ Plots the trajectory of the slithering sphericalJoint
    as simulation progresses
    """

    def __init__(
        self, input_folder, output_folder, save_file, force_flag, display_flag
    ):
        super(HingeJointPlotter, self).__init__(
            input_folder, output_folder, save_file, force_flag, display_flag
        )

        # Overrides changes
        self.ax.set_xlim(-200, 200)
        self.ax.set_ylim(-150, 250)
        self.ax.set_zlim(00, 400)
        set_axes_equal(self.ax)


class FixedJointPlotter(SphericalJointPlotter):
    """ Plots the trajectory of the slithering sphericalJoint
    as simulation progresses
    """

    def __init__(
        self, input_folder, output_folder, save_file, force_flag, display_flag
    ):
        super(FixedJointPlotter, self).__init__(
            input_folder, output_folder, save_file, force_flag, display_flag
        )

        # # Overrides changes
        # self.ax.set_xlim(-50, 50)
        # self.ax.set_ylim(-150, 250)
        # self.ax.set_zlim(300, 400)
        # set_axes_equal(self.ax)


class PullingMusclePlotter(ThreeDimensionalPlotter):
    """ Plots the trajectory of the slithering sphericalJoint
    as simulation progresses
    """

    def __init__(
        self, input_folder, output_folder, save_file, force_flag, display_flag
    ):
        # Body and wireframe color
        self.file_metadata = [
            # ("rod1", ".txt", (31 / 255, 119 / 255, 180 / 255)),
            ("rod1", ".txt", to_rgb("xkcd:bluish")),
            # ("rod2", ".txt", (200 / 255, 0 / 255, 0 / 255)),
            ("rod2", ".txt", to_rgb("xkcd:reddish")),
            # ("rod3", ".txt", (50 / 255, 200 / 255, 80 / 255)),
            ("rod3", ".txt", to_rgb("xkcd:greenish")),
        ]

        super(PullingMusclePlotter, self).__init__(
            input_folder, output_folder, save_file, force_flag, display_flag
        )

        # Make any other changes to the figure here
        # Data-dependent maybe
        self.ax.set_xlim(-50, 50)
        self.ax.set_ylim(-100, 200)
        self.ax.set_zlim(50, 350)
        # self.ax.set_ylim(-0.125, 0.125)
        # self.ax.set_zlim(-0.2, 0.2)
        set_axes_equal(self.ax)


class SnakePlotter(ThreeDimensionalPlotter):
    """ Plots the trajectory of the slithering snake
    as simulation progresses
    """

    def __init__(
        self, input_folder, output_folder, save_file, force_flag, display_flag
    ):
        # Body and wireframe color
        self.file_metadata = [("rod1", ".txt", to_rgb("xkcd:bluish"))]

        super(SnakePlotter, self).__init__(
            input_folder, output_folder, save_file, force_flag, display_flag
        )

        # Make any other changes to the figure here
        # Data-dependent maybe
        self.ax.set_xlim(-3, 0.8)
        self.ax.set_ylim(-0.25, 0.25)
        self.ax.set_zlim(-0.05, 0.05)
        set_axes_equal(self.ax)

    def plot(self, ax, data, color=(31 / 255, 119 / 255, 180 / 255), i_seq=None):
        super(SnakePlotter, self).plot(ax, data, color, i_seq)

        # After the last possible dynamic update plot
        if i_seq == 0:  # ie. i_seq == len(self.metadata)-1
            if len(ax.collections) > 2:  # 2 * i_seq + 2
                # surfae alrady plotted, why plot it again?
                pass
            else:
                # * unpakcs arguments
                x = np.linspace(*ax.get_xlim(), 7)
                y = np.linspace(*ax.get_ylim(), 7)
                X, Y = np.meshgrid(x, y)
                Z = 0.0 * X
                ax.grid(False)
                ax.plot_wireframe(
                    X, Y, Z, rstride=1, cstride=1, color="darkgrey", linewidth=1
                )

    # def plot(self, ax, data, color=(31 / 255, 119 / 255, 180 / 255)):
    #     """ Plots snake data
    #     """
    #     super(SnakePlotter, self).plot(ax, data, color)

    #     # self.ax.set_xlim(0, 0.125)
    #     # self.ax.set_ylim(-0.125, 0.125)
    #     # self.ax.set_zlim(-0.2, 0.2)
    #     set_axes_equal(self.ax)


class HelicalBucklingPlotter(ThreeDimensionalPlotter):
    """ Plots the trajectory of the slithering snake
    as simulation progresses
    """

    def __init__(
        self, input_folder, output_folder, save_file, force_flag, display_flag
    ):
        # Body and wireframe color
        self.file_metadata = [("rod1", ".txt", to_rgb("xkcd:bluish"))]

        super(HelicalBucklingPlotter, self).__init__(
            input_folder, output_folder, save_file, force_flag, display_flag
        )

        # Make any other changes to the figure here
        # Data-dependent maybe
        # self.ax.set_xlim(-3, 0.8)
        # self.ax.set_ylim(-0.25, 0.25)
        self.ax.set_zlim(-50.0, 50.0)
        set_axes_equal(self.ax)

    def plot(self, ax, data, color=(31 / 255, 119 / 255, 180 / 255), i_seq=None):
        super(HelicalBucklingPlotter, self).plot(ax, data, color, i_seq)


class TwoDimensionalPlotter:
    """ Plot class for classical 2D line plots"""

    def __init__(
        self, input_folder, output_folder, save_file, force_flag, display_flag
    ):

        """ Figure attributes """
        self.fig = plt.figure(figsize=FigProperties.figsize())

        if not hasattr(self, "file_metadata"):
            #     print("hi")
            # if self.file_metadata is None:
            self.file_metadata = [("Flagella", ".txt", to_rgb("xkcd:bluish"))]

        # Also total number of suplots
        ntypes_files = len(self.file_metadata)
        # Number of rows dependent
        n_rows = 1 if ntypes_files < 4 else 2
        # Number of columns,
        n_columns = ntypes_files // n_rows
        n_columns += ntypes_files % n_rows

        # Create a Position index
        pindex = range(1, ntypes_files + 1)

        # Create a bunch of axes
        self.axes = [self.fig.add_subplot(n_rows, n_columns, k) for k in pindex]

        # Set some common attributes
        for axis in self.axes:
            axis.grid(True)
            # ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

        """ Save file attributes """
        self.output_folder = output_folder
        self.savefile_name, self.savefile_ext = os.path.splitext(save_file)

        # Assume its a pdf if no extension given by default
        if not self.savefile_ext:
            self.savefile_ext = ".pdf"

        """ Source file attributes and search """

        # Guaranteed to be sorted
        self.files_to_be_processed = [None for i in range(ntypes_files)]
        for index, (prefix, suffix, _) in enumerate(self.file_metadata):
            temp, _ = detect_files(input_folder, prefix, suffix)
            self.files_to_be_processed[index] = temp.pop()

        self.force_flag = force_flag
        self.display_flag = display_flag

        print("{} initialized".format(type(self).__name__))

    def process(self):
        """ Loads all data, plots them and stores them into appropriately
        named figures.
        """

        # Load color metadata first as list of rgb tuples
        colors = [color for (_, _, color) in self.file_metadata]

        if self.display_flag:
            # Turn on interactive mode to persist figure
            plt.ion()
            # Show figure after persist
            plt.show()

        print("Starting processing, press Ctrl+C to quit at any time")
        print("(Maybe more than once)")
        while True:
            try:
                for i_ax, (src_file, axis, color) in enumerate(
                    zip(self.files_to_be_processed, self.axes, colors)
                ):
                    # Defaults loads to (ndata, 4) rather than (4,ndata)
                    data = np.loadtxt(src_file).T
                    self.plot(axis, data, color, i_ax)

                    # Autoscaling important as we dynamically update now
                    axis.relim()
                    axis.autoscale_view(True, True, True)

                self.fig.canvas.draw()

                filename = "{name}{ext}".format(
                    name=self.savefile_name, ext=self.savefile_ext
                )

                self.fig.savefig(
                    os.path.join(self.output_folder, filename), dpi=FigProperties.dpi
                )

                if self.display_flag:
                    # Bigger pause to prevent too much looping
                    plt.pause(2.0)
                else:
                    pass

            except KeyboardInterrupt:
                print("Finished processing")
                raise IOError

        # if self.display_flag:
        #     # Show only the last drawn figure to the user
        #     plt.show()

    def plot(self, axis, data, color=(31 / 255, 119 / 255, 180 / 255), i_ax=None):
        """ Plots snake data
        """
        if axis.lines:
            for line in axis.lines:
                line.set_data(data[0, :], data[1, :])
        else:
            axis.plot(data[0, :], data[1, :], color=color, marker="o", linewidth=2)
        # axis.plot(data[0, :], data[1, :], color=color, marker="o", linewidth=2)

    def animate(self):
        """ Animates using ffmpeg the figures created by savefig
        """
        pass


class TimoshenkoPlotter(TwoDimensionalPlotter):
    """Plots flagella quantities"""

    # pylint : disable=too-many-arguments
    def __init__(
        self, input_folder, output_folder, save_file, force_flag, display_flag
    ):
        # pylint : enable=too-many-arguments

        self.file_metadata = [("timoshenko_final_shape", ".txt", to_rgb("xkcd:bluish"))]

        super(TimoshenkoPlotter, self).__init__(
            input_folder, output_folder, save_file, force_flag, display_flag
        )

        # Setting labels here
        self.axes[0].set_xlabel("X (m)")
        self.axes[0].set_ylabel("Y (m)")
        self.axes[0].set_title("Centerline deflection")

        print("TimoshenkoPlotter initialized!")

    def plot(self, ax, data, color=(31 / 255, 119 / 255, 180 / 255), i_ax=None):
        """ Plots timoshenko data
        """
        # Timoshenko data needs correction with the first instant
        ax.plot(
            data[0, :],
            data[1, :],
            color=color,
            marker="o",
            linewidth=2,
            label="simulation",
        )
        ax.plot(
            data[3, :],
            data[4, :],
            color="k",
            linestyle="dashed",
            linewidth=2,
            label="analytical",
        )
        ax.legend()

        # In this case as soon as we plot, we can exit
        plt.pause(10.0)
        import sys

        sys.exit(0)


class ElbowPlotter(TwoDimensionalPlotter):
    """ Plots the muscle force output of the elbow
    as simulation progresses
    """

    # pylint : disable=too-many-arguments
    def __init__(
        self, input_folder, output_folder, save_file, force_flag, display_flag
    ):
        # pylint : enable=too-many-arguments

        self.file_metadata = [("velocity", ".txt", to_rgb("xkcd:bluish"))]

        super(ElbowPlotter, self).__init__(
            input_folder, output_folder, save_file, force_flag, display_flag
        )

        # Make any other changes to the figure here
        # Data-dependent maybe

        # Setting labels here
        self.axes[0].set_xlabel("T (s)")
        self.axes[0].set_ylabel("Elbow angle (degrees)")
        self.axes[0].set_title("Angle vs T")

        # self.ax.set_xlim(0, 1)
        # self.ax.set_ylim(0, 1)

    def plot(self, ax, data, color=(31 / 255, 119 / 255, 180 / 255), i_ax=None):
        """ Plots elbow data
        """
        # Flagella data needs correction with the first instant
        time = data[0, :]
        angle = data[1, :] * 180.0 / np.pi
        # print(time, pos)

        if ax.lines:
            for line in ax.lines:
                line.set_data(time, angle)
        else:
            ax.plot(time, angle, color=color, marker="o", linewidth=2)


# class SnakeVelocityPlotter(TwoDimensionalPlotter):
#     """ Plots the muscle force output of the elbow
#     as simulation progresses
#     """

#     # pylint : disable=too-many-arguments
#     def __init__(
#         self, input_folder, output_folder, save_file, force_flag, display_flag
#     ):
#         # pylint : enable=too-many-arguments

#         # self.file_metadata = [("velocity", ".txt", (31 / 255, 119 / 255, 180 / 255))]

#         super(SnakeVelocityPlotter, self).__init__(
#             input_folder, output_folder, save_file, force_flag, display_flag
#         )

#         # Make any other changes to the figure here
#         # Data-dependent maybe

#         # Setting labels here
#         self.axes[0].set_xlabel("T (s)")
#         self.axes[0].set_ylabel("Velocity (m/s)")
#         self.axes[0].set_title("Velocity vs T")

#         # self.ax.set_xlim(0, 1)
#         # self.ax.set_ylim(0, 1)

#     def plot(self, ax, data, color=(31 / 255, 119 / 255, 180 / 255), i_ax=None):
#         """ Plots muscularsnake data
#         """
#         # Muscularsnake data plotted from the second instant
#         # The first instant is just the initialization which is
#         # discarded to give time for the structure to relax (ie contact etc.)

#         # Muscularsnake data plotted with reference to the second instant
#         # The first is initial time
#         time = data[0, :]
#         x_velocity = data[3, :]
#         y_velocity = data[4, :]

#         velocity_dir = np.array([-2.385, -0.089])
#         velocity_dir /= norm(velocity_dir)

#         fwd_velocity = x_velocity * velocity_dir[0] + y_velocity * velocity_dir[1]

#         binorm = np.array([0.089, -2.385])
#         binorm /= norm(binorm)

#         lat_velocity = x_velocity * binorm[0] + y_velocity * binorm[1]
#         mag_vel_new = fwd_velocity**2 + lat_velocity**2
#         mag_vel_orig = x_velocity**2 + y_velocity**2

#         if ax.lines:
#             ax.lines[0].set_data(time, -fwd_velocity)
#             ax.lines[1].set_data(time, lat_velocity)
#             # ax.lines[0].set_data(time, mag_vel_orig)
#             # ax.lines[1].set_data(time, mag_vel_new)
#         else:
#             ax.plot(
#                 time,
#                 fwd_velocity,
#                 color=color,
#                 marker="o",
#                 linewidth=2,
#                 label="forward",
#             )
#             ax.plot(
#                 time,
#                 lat_velocity,
#                 # color="xkcd : deep red",
#                 # c="xkcd: deep red",
#                 c = 'k',
#                 marker="o",
#                 linewidth=2,
#                 label="lateral",
#             )
#             ax.legend()


class FlagellaPlotter(TwoDimensionalPlotter):
    """Plots flagella quantities"""

    # pylint : disable=too-many-arguments
    def __init__(
        self, input_folder, output_folder, save_file, force_flag, display_flag
    ):
        # pylint : enable=too-many-arguments

        # self.file_metadata = [("Flagella", ".txt", (31 / 255, 119 / 255, 180 / 255))]

        super(FlagellaPlotter, self).__init__(
            input_folder, output_folder, save_file, force_flag, display_flag
        )

        # Make any other changes to the figure here
        # Data-dependent maybe

        # Setting labels here
        self.axes[0].set_xlabel("Time (s)")
        self.axes[0].set_ylabel("Position (micrometer)")
        self.axes[0].set_title("x CoM position vs Time")

    def plot(self, ax, data, color=(31 / 255, 119 / 255, 180 / 255), i_ax=None):
        """ Plots flagella data
        """
        # Flagella data needs correction with the first instant
        time = data[0, :]
        pos = data[1, :] - data[1, 0]
        # print(time, pos)

        if ax.lines:
            for line in ax.lines:
                line.set_data(time, -1000.0 * pos)
        else:
            ax.plot(time, -1000 * pos, color=color, marker="o", linewidth=2)


class WalkerPlotter(TwoDimensionalPlotter):
    """Plots flagella quantities"""

    # pylint : disable=too-many-arguments
    def __init__(
        self, input_folder, output_folder, save_file, force_flag, display_flag
    ):
        # pylint : enable=too-many-arguments

        # self.file_metadata = [("Flagella", ".txt", (31 / 255, 119 / 255, 180 / 255))]

        super(WalkerPlotter, self).__init__(
            input_folder, output_folder, save_file, force_flag, display_flag
        )

        # Make any other changes to the figure here
        # Data-dependent maybe

        # Setting labels here
        self.axes[0].set_xlabel("Time (s)")
        self.axes[0].set_ylabel("Displacement (millimeter)")
        self.axes[0].set_title("Displacement vs Time")

    def plot(self, ax, data, color=(31 / 255, 119 / 255, 180 / 255), i_ax=None):
        """ Plots walker data
        """
        # Walker data plotted from the second instant
        # The first instant is just the initialization which is
        # discarded to give time for the structure to relax (ie contact etc.)

        # Walker data plotted with reference to the second instant
        # The first is initial time
        time = data[0, 1:]
        pos = data[1, 1:] - data[1, 1]

        if ax.lines:
            for line in ax.lines:
                line.set_data(time, -pos)
        else:
            # - in pos to correct for direction
            ax.plot(time, -pos, color=color, marker="o", linewidth=2)


class SnakeVelocityPlotter(TwoDimensionalPlotter):
    """Plots snake quantities"""

    # pylint : disable=too-many-arguments
    def __init__(
        self, input_folder, output_folder, save_file, force_flag, display_flag
    ):
        # pylint : enable=too-many-arguments

        # self.file_metadata = [("Flagella", ".txt", (31 / 255, 119 / 255, 180 / 255))]

        super(SnakeVelocityPlotter, self).__init__(
            input_folder, output_folder, save_file, force_flag, display_flag
        )

        # Setting labels here
        self.axes[0].set_xlabel("T")
        self.axes[0].set_ylabel("V (m/s)")
        self.axes[0].set_title("Velocity vs T")

    def plot(self, ax, data, color=(31 / 255, 119 / 255, 180 / 255), i_ax=None):
        """ Plots muscularsnake data
        """
        # Muscularsnake data plotted from the second instant
        # The first instant is just the initialization which is
        # discarded to give time for the structure to relax (ie contact etc.)

        # Muscularsnake data plotted with reference to the second instant
        # The first is initial time
        time = data[0, :]
        fwd_velocity = data[3, :]
        lat_velocity = data[4, :]

        if ax.lines:
            ax.lines[0].set_data(time, fwd_velocity)
            ax.lines[1].set_data(time, lat_velocity)
        else:
            ax.plot(
                time,
                fwd_velocity,
                color=color,
                marker="o",
                linewidth=2,
                label="forward",
            )
            ax.plot(
                time,
                lat_velocity,
                c="xkcd:reddish",
                marker="o",
                linewidth=2,
                label="lateral",
            )
            ax.legend()

        self.file_metadata = [("timoshenko_final_shape", ".txt", to_rgb("xkcd:bluish"))]

        super(TimoshenkoPlotter, self).__init__(
            input_folder, output_folder, save_file, force_flag, display_flag
        )

        # Setting labels here
        self.axes[0].set_xlabel("X (m)")
        self.axes[0].set_ylabel("Y (m)")
        self.axes[0].set_title("Centerline deflection")

        print("TimoshenkoPlotter initialized!")

    def plot(self, ax, data, color=(31 / 255, 119 / 255, 180 / 255), i_ax=None):
        """ Plots timoshenko data
        """
        # Timoshenko data needs correction with the first instant
        ax.plot(
            data[0, :],
            data[1, :],
            color=color,
            marker="o",
            linewidth=2,
            label="simulation",
        )
        ax.plot(
            data[3, :],
            data[4, :],
            color="k",
            linestyle="dashed",
            linewidth=2,
            label="analytical",
        )
        ax.legend()

        # In this case as soon as we plot, we can exit
        plt.pause(10.0)
        import sys

        sys.exit(0)


class HelicalPhiPlotter(TwoDimensionalPlotter):
    """Plots helical phi """

    # pylint : disable=too-many-arguments
    def __init__(
        self, input_folder, output_folder, save_file, force_flag, display_flag
    ):
        # pylint : enable=too-many-arguments

        self.file_metadata = [("helix_0100_shape", ".txt", to_rgb("xkcd:bluish"))]

        super(HelicalPhiPlotter, self).__init__(
            input_folder, output_folder, save_file, force_flag, display_flag
        )
        # length of rod
        self.L = 100.0

        # Setting labels here
        self.axes[0].set_xlabel("s - L/2")
        self.axes[0].set_ylabel("phi")
        self.axes[0].set_title("phi vs s")

    def envelope(self, arg_pos):
        """
        Given points, computes the arc length and envelope of the curve
        """
        n_points = arg_pos.shape[1]

        # Computes the direction in which the rod points
        # in our cases it should be the z-axis
        rod_direction = arg_pos[:, -1] - arg_pos[:, 0]
        rod_direction /= norm(rod_direction, ord=2, axis=0)

        # Compute local tangent directions
        tangent_s = np.diff(arg_pos, n=1, axis=-1)  # x_(i+1)-x(i)
        length_s = norm(tangent_s, ord=2, axis=0)
        tangent_s /= length_s

        # Dot product with direction is cos_phi, see RSOS
        cos_phi_s = np.einsum("ij,i->j", tangent_s, rod_direction)

        # Compute phi-max now
        phi = np.arccos(cos_phi_s)
        cos_phi_max = np.cos(np.max(phi))

        # Return envelope and arclength
        envelope = (cos_phi_s - cos_phi_max) / (1.0 - cos_phi_max)
        # -0.5 * length accounts for the element/node business
        arclength = np.cumsum(length_s) - 0.5 * length_s[0]

        return arclength, envelope

    def analytical_solution(self):
        """ Gives the analytical solution of the helicalbuckling case
        """
        # Physical parameters, set from the simulation
        B = 1.345
        C = 0.789
        gamma = C / B
        R = 27.0 * 2.0 * np.pi
        d = 0.03
        D = d * self.L
        nu = 1.0 / gamma - 1.0

        # These are magic constants, but you can obtain them by solving
        # this equation (accoring to matlab syntax)
        # syms x y
        # S = vpasolve([d == sqrt(16/y*(1-x*x/(4*y))), R == x/gamma+4*acos(x/(2*sqrt(y)))], [x, y]);
        # moment = double(S.x); # dimensionless end moment
        # tension = double(S.y); # dimensionless end torque
        # This comes from  Eqs. 14-15 of "Writhing instabilities of twisted rods: from
        # infinite to finite length", 2001
        # We did not want to introduce sympy dependency here, so we decided to hardcode
        # the solutions instead
        moment = 98.541496171190744
        tension = 2.900993205792131e3

        # Compute maximum envelope angle according to Eq. 13 of "Writhing
        # instabilities of twisted rods: from infinite to finite length", 2001
        thetaMax = np.arccos(moment * moment / (2.0 * tension) - 1.0)

        # Compute actual end torque and tension according to "Writhing
        # instabilities of twisted rods: from infinite to finite length", 2001
        M = moment * B / self.L
        T = tension * B / (self.L * self.L)

        # Compute dimensionless load according to Eq. 30 of "Helical and localised
        # buckling in twisted rods: a unified analysis of the symmetric case", 2000
        m = M / np.sqrt(B * T)

        # Setup for analytical curve calculation
        s = np.linspace(-0.5, 0.5, 10000)
        t = T * self.L * self.L / (4 * np.pi * np.pi * B)
        mz = M * self.L / (2 * np.pi * B)
        root = np.sqrt(4 * t - mz * mz)

        # This is the analytical curve computed
        # according to Eqs. 27 and 52 of
        # "Instability and self-contact phenomena in the writhing of clamped rods",
        # 2003
        xs = (
            1.0
            / (2.0 * np.pi * t)
            * root
            * np.sin(mz * np.pi * s)
            / np.cosh(np.pi * s * root)
        )
        ys = (
            -1.0
            / (2.0 * np.pi * t)
            * root
            * np.cos(mz * np.pi * s)
            / np.cosh(np.pi * s * root)
        )
        zs = s - 1.0 / (2.0 * np.pi * t) * root * np.tanh(np.pi * s * root)
        pos = np.vstack((xs, ys, zs)) * self.L
        return self.envelope(pos)

    def plot(self, ax, data, color=(31 / 255, 119 / 255, 180 / 255), i_ax=None):
        """ Plots muscularsnake data
        """
        # Muscularsnake data plotted from the second instant
        # The first instant is just the initialization which is
        # discarded to give time for the structure to relax (ie contact etc.)

        # Muscularsnake data plotted with reference to the second instant
        # The first is initial time
        if ax.lines:
            pass
        else:
            analytical_centerline, analytical_envelope = self.analytical_solution()
            num_centerline, num_envelope = self.envelope(data)

            ax.plot(
                num_centerline - 0.5 * self.L,
                num_envelope,
                color=color,
                marker="o",
                linewidth=2,
                label="numerical",
            )
            ax.plot(
                analytical_centerline - 0.5 * self.L,
                analytical_envelope,
                c="black",
                linestyle="--",
                linewidth=1,
                label="analytical",
            )
            ax.legend()


class WingPlotter(TwoDimensionalPlotter):
    """Plots wing quantities"""

    # pylint : disable=too-many-arguments
    def __init__(
        self, input_folder, output_folder, save_file, force_flag, display_flag
    ):
        # pylint : enable=too-many-arguments

        # self.file_metadata = [("Flagella", ".txt", (31 / 255, 119 / 255, 180 / 255))]

        super(WingPlotter, self).__init__(
            input_folder, output_folder, save_file, force_flag, display_flag
        )

        # Setting labels here
        self.axes[0].set_xlabel("T (s)")
        self.axes[0].set_ylabel("Angle (degrees)")
        self.axes[0].set_title("Angle vs T")

    def plot(self, ax, data, color=(31 / 255, 119 / 255, 180 / 255), i_ax=None):
        """ Plots wing data
        """
        # Wing data plotted from the second instant
        # The first instant is just the initialization which is
        # discarded to give time for the structure to relax (ie contact etc.)

        # Wing data plotted with reference to the second instant
        # The first is initial time
        time = (data[0, 1:] - 0.125) / (
            0.38
        )  # Accounts for initialization time and period
        dv_angle = data[1, 1:]  # dorsoventral
        ap_angle = data[2, 1:]  # anterio-posterior
        elb_angle = data[3, 1:]  # elbow-angle

        if ax.lines:
            ax.lines[0].set_data(time, dv_angle)
            ax.lines[1].set_data(time, ap_angle)
            ax.lines[2].set_data(time, elb_angle)
        else:
            ax.plot(
                time,
                dv_angle,
                color=color,
                marker="o",
                linewidth=2,
                label="dorsoventral",
            )
            ax.plot(
                time,
                ap_angle,
                color="xkcd:reddish",
                marker="o",
                linewidth=2,
                label="anterio-posterior",
            )
            ax.plot(
                time,
                elb_angle,
                color="xkcd:greenish",
                marker="o",
                linewidth=2,
                label="elbow",
            )
            ax.legend()


def parse_args(args=None):
    """Parse arguments from commandline"""
    parser = argparse.ArgumentParser(
        description=sys.modules[__name__].__doc__, formatter_class=CustomFormatter
    )

    parser.add_argument(
        "-c",
        "--case",
        choices=[
            "timoshenkobeam",
            "helicalbuckling",
            "sphericaljoint",
            "hingejoint",
            "fixedjoint",
            "pullingmuscle",
            "snake",
            "elbow",
            "flagella",
            "walker",
            "muscularsnake",
            "wing",
        ],
        help="simulation case whose output needs to be seen",
        type=str,
    )

    parser.add_argument(
        "-o",
        "--output",
        metavar="OUTFILE",
        help="if enabled, stores the images as OUTFILE_{1,2,...,N}.pdf,\
        in the scripts/pyprocessed_<case> directory",
        default="out",
        type=str,
    )

    parser.add_argument(
        "-p",
        "--path",
        metavar="OUTPATH",
        help="path to store output files, is default created to\
        pyprocessed_case if not initialized",
        default="./pyprocessed",
        type=str,
    )

    parser.add_argument(
        "-f",
        "--force",
        help="force rewrite of any previously saved images",
        action="store_true",
    )

    parser.add_argument(
        "-a",
        "--animate",
        help="force collate pngs (ffmpeg) after saving them",
        action="store_true",
    )

    parser.add_argument(
        "-n",
        "--nodisp",
        help="do not render images on the screen but only save (faster!)",
        action="store_true",
    )

    if len(args) < 2:
        parser.print_help()
        sys.exit(42)

    return parser.parse_args(args)


def find_results(folder_name):
    """ Given a folder name, traverse the current dir to find the folder,if not
    found travel its parent, and so on until root.
    Credits: https://stackoverflow.com/a/37560251
    """
    cur_dir = os.getcwd()

    while True:
        parent_dir = os.path.dirname(cur_dir)
        test_dir = os.path.join(cur_dir, folder_name)
        #
        if os.path.isdir(test_dir):  # pylint : disable=no-else-return
            return test_dir
        else:
            if cur_dir == parent_dir:  # if dir is root dir
                raise FileNotFoundError("Folder {} not found".format(folder_name))
            else:
                cur_dir = parent_dir


def main(argv):
    """ main function coordinating output
    """

    """Parse opts"""
    options = parse_args(argv)

    case_name = options.case
    output_file = options.output
    output_folder = options.path
    force_flag = options.force
    animate_flag = options.animate
    display_flag = not options.nodisp

    """Process flags"""
    input_folder = find_results("run_" + case_name)

    output_folder = os.path.abspath(os.path.join(output_folder, case_name))
    # Lazy to implement version safe code
    os.makedirs(output_folder, exist_ok=True)

    tre_d_plotter = None
    two_d_plotter = None

    two_display_flag = display_flag
    three_display_flag = display_flag

    """Decide and process"""
    if case_name == "timoshenkobeam":
        two_d_plotter = TimoshenkoPlotter
    elif case_name == "helicalbuckling":
        two_d_plotter = HelicalPhiPlotter
        tre_d_plotter = HelicalBucklingPlotter
    elif case_name == "elbow":
        two_d_plotter = ElbowPlotter
    elif case_name == "flagella":
        two_d_plotter = FlagellaPlotter
    elif case_name == "walker":
        two_d_plotter = WalkerPlotter
    elif case_name == "muscularsnake":
        two_d_plotter = SnakeVelocityPlotter
    elif case_name == "wing":
        two_d_plotter = WingPlotter
    elif case_name == "sphericaljoint":
        tre_d_plotter = SphericalJointPlotter
    elif case_name == "hingejoint":
        tre_d_plotter = HingeJointPlotter
    elif case_name == "fixedjoint":
        tre_d_plotter = FixedJointPlotter
    elif case_name == "pullingmuscle":
        tre_d_plotter = PullingMusclePlotter
    elif case_name == "snake":
        two_d_plotter = SnakeVelocityPlotter
        tre_d_plotter = SnakePlotter

    if tre_d_plotter is None:
        tre_d_plotter = DummyPlotter
        three_display_flag = False
    if two_d_plotter is None:
        two_d_plotter = DummyPlotter
        two_display_flag = False

    # Order important here, people are interested in seeing the plots first
    # to see whether something happened or not
    plotters = [
        tre_d_plotter(
            input_folder, output_folder, output_file, force_flag, three_display_flag
        ),
        two_d_plotter(
            input_folder, output_folder, output_file, force_flag, two_display_flag
        ),
    ]

    # pylint: disable=expression-not-assigned
    if not animate_flag:
        [plotter.process() for plotter in plotters]
    elif force_flag and animate_flag:
        [plotter.process() for plotter in plotters]
        [plotter.animate() for plotter in plotters]
    else:
        # Assume images have been done at output folder, you can just animate them now
        [plotter.animate() for plotter in plotters]
    # pylint: enable=expression-not-assigned


if __name__ == "__main__":
    main(sys.argv[1:])
