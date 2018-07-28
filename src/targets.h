// Licensed under a 3-clause BSD style license - see LICENSE.rst

#ifndef TARGETS_H
#define TARGETS_H

#include <cstdint>

#include <string>
#include <vector>
#include <map>
#include <set>
#include <memory>
#include <algorithm>
#include <exception>
#include <sstream>

#include <htmTree.h>
#include <kdTree.h>

#include <utils.h>
#include <tiles.h>


namespace fiberassign {

// Allowed target types

#define TARGET_TYPE_SCIENCE 1
#define TARGET_TYPE_STANDARD 2
#define TARGET_TYPE_SKY 4
#define TARGET_TYPE_SAFE 8

std::string target_string(uint8_t type);


// This simple class represents the properties of a single target.
// This is only used internally and is not exposed to Python.

class Target {

    public :

        typedef std::shared_ptr <Target> pshr;

        Target();

        Target(
            int64_t tid,
            double tra,
            double tdec,
            int32_t tobs_remain,
            int32_t tpriority,
            double tsubpriority,
            int32_t tobscond,
            uint8_t ttype
        );

        int64_t id;
        double ra;
        double dec;
        int32_t obs_remain;
        int32_t priority;
        double subpriority;
        int32_t obscond;
        uint8_t type;

        bool is_science() const;
        bool is_standard() const;
        bool is_sky() const;
        bool is_safe() const;
        bool is_type(uint8_t t) const;

};


// This class holds the information for multiple targets

class Targets : public std::enable_shared_from_this <Targets> {

    public :

        typedef std::shared_ptr <Targets> pshr;

        Targets();

        void append (
            std::vector <int64_t> const & id,
            std::vector <double> const & ra,
            std::vector <double> const & dec,
            std::vector <int32_t> const & obs_remain,
            std::vector <int32_t> const & priority,
            std::vector <double> const & subpriority,
            std::vector <int32_t> const & obscond,
            std::vector <uint8_t> const & type
        );

        std::map <int64_t, Target> data;
        std::set <int32_t> science_classes;

};


// This is an HTM tree that indexes into an existing list of objects.

typedef struct {
    int64_t id;
    double nhat[3];
} TreePoint;

class TargetTree : public std::enable_shared_from_this <TargetTree> {

    public :

        typedef std::shared_ptr <TargetTree> pshr;

        TargetTree(Targets::pshr objs, double min_tree_size = 0.01);

        void near(double ra_deg, double dec_deg, double radius_rad,
            std::vector <int64_t> & result) const;

    private :

        // Helper vector with just the data we need for the tree.
        std::vector <TreePoint> treelist_;

        // Wrap in a unique pointer so that resources are freed automatically
        // on destruction.
        std::unique_ptr < htmTree <TreePoint> > tree_;

        // min tree size parameter.
        double mintreesz_;

};


// Data for one point of the KD tree.

typedef struct {
    int64_t id;
    double pos[2];
} KdTreePoint;

// Class holding the object IDs available for each tile and fiber.

class TargetsAvailable : public std::enable_shared_from_this <TargetsAvailable> {

    public :

        typedef std::shared_ptr <TargetsAvailable> pshr;

        TargetsAvailable(Targets::pshr objs, Tiles::pshr tiles,
                         TargetTree::pshr tree);

        Tiles::pshr tiles() const;

        std::map <int32_t, std::vector <int64_t> > tile_data(int32_t tile) const;

        std::map <int32_t, std::map <int32_t, std::vector <int64_t> > > data;

    private :

        Tiles::pshr tiles_;

};


// Class holding the tile / fibers available for each target ID.

class FibersAvailable : public std::enable_shared_from_this <FibersAvailable> {

    public :

        typedef std::shared_ptr <FibersAvailable> pshr;

        FibersAvailable(TargetsAvailable::pshr tgsavail);

        std::vector <std::pair <int32_t, int32_t> >
            target_data(int64_t target) const;

        std::map < int64_t, std::vector < std::pair <int32_t, int32_t> > > data;

};





}
#endif
