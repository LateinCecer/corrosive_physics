# Corrosive Physics

Corrosive Physics is an open source real-time physics engine written in Rust. It
aims to provide cutting-edge performance, especially for rigid body dynamics.
Support for soft-bodies is not yet planned, but may be included in the future.


# Data Structures

The data structures implemented so far are a basic BVH and a TLAS. These may be
used for collision queries.

# Rendering Frontend

For testing purposes, having a rendering frontend would be desirable. Since this
crate is written in Rust, it makes sense to use a Rust-based solution for this.
As such, integration of the Corrosive Physics engine with 
[Bevys](https://crates.io/crates/bevy) ECS is planned.

# License

This software is distributed under the MIT license (see LICENSE file in the root
folder of this project).

# Contributions

If you like this project and want to contribute, feel free to create a pull request
:)
