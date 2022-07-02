use bevy::prelude::*;
use nalgebra::Vector3;
use corrosive_physics::engine::{PhysicsEngine};
use corrosive_physics::system::object::{PhyEntity, PhyEntityID};
use corrosive_physics::volume::BVIntersector;
use corrosive_physics::volume::tlas::TLASElement;

fn main() {
    println!("String test case 'Cubes'...");
    App::new()
        .insert_resource(Msaa { samples: 4 })
        .add_plugins(DefaultPlugins)
        .add_startup_system(setup)
        .add_system(update)
        .run();
}

#[derive(Component)]
struct Rotator;

fn update(
    time: Res<Time>,
    mut query: Query<(&PhyEntityID, &mut Transform)>
) {
    let mut engine = PhysicsEngine::global_mut();

    for (id, mut trans) in query.iter_mut() {
        // eprintln!("ids {:?}  @  {:?}", id.entity_id, trans.translation);

        if id.entity_id != 0 {
            // query for potential colliders
            let colliders = engine.query_colliders(id.clone());
            // let floor = &engine[PhyEntityID { world_id: 0, chunk_id: 0, entity_id: 0 }];

            if colliders.is_empty() || (colliders.len() == 1 && colliders[0].id.entity_id == id.entity_id) {
                // if !engine.world.nodes()[1].aabb().intersects(engine[id.clone()].bounding_volume()) {
                // update
                let entity: &mut PhyEntity<f64> = &mut engine[id.clone()];
                entity.is.integrate(time.delta_seconds_f64());
                entity.sync();

                // refit TLAS to the updated bounds (faster than a full rebuild)
                engine.world.refit();
            }
        }



        // sync
        let entity: &PhyEntity<f64> = &engine[id.clone()];
        let transform: &mut Transform = &mut trans;

        transform.translation.x = entity.is.state.pos.x as f32;
        transform.translation.y = entity.is.state.pos.y as f32;
        transform.translation.z = entity.is.state.pos.z as f32;

        transform.rotation = Quat::from_xyzw(
            entity.is.state.rot.i as f32,
            entity.is.state.rot.j as f32,
            entity.is.state.rot.k as f32,
            entity.is.state.rot.w as f32,
        );

        transform.scale.x = entity.is.state.scale.x as f32;
        transform.scale.y = entity.is.state.scale.y as f32;
        transform.scale.z = entity.is.state.scale.z as f32;
    }

    // rebuild the tree properly for the next tick
    engine.world.build();
}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    let cube_handle = meshes.add(Mesh::from(shape::Cube { size: 1.0 }));
    let cube_material_handle = materials.add(
        StandardMaterial {
            base_color: Color::rgb(0.8, 0.7, 0.6),
            .. default()
        }
    );

    let floor_handle = meshes.add(Mesh::from(shape::Box {
        min_x: -5.0,
        max_x: 5.0,
        min_y: -0.5,
        max_y: 0.5,
        min_z: -5.0,
        max_z: 5.0
    }));
    let floor_material_handle = materials.add(
        StandardMaterial {
            base_color: Color::rgb(0.3, 0.5, 0.8),
            .. default()
        }
    );




    // create engine and physical shadows of renderable objects
    let mut engine = PhysicsEngine::<f64>::new();
    let mut count = 0usize;
    let floor_id = PhyEntityID {
        world_id: 0,
        chunk_id: 0,
        entity_id: count,
    };
    count += 1;

    let mut floor = PhyEntity::cube(
        floor_id.clone(), Vector3::new(20.0, 1.0, 20.0));

    floor.is.state.pos = Vector3::new(0.0, 0.0, 0.0);
    floor.is.momentum = Vector3::new(0.0, 0.0, 0.0);
    floor.sync();
    engine.world.blas_mut().push(floor);


    let spacing = 2.0;
    for y in 0..10 {
        for x in 0..6 {
            for z in 0..6 {

                let cube_id = PhyEntityID {
                    world_id: 0,
                    chunk_id: 0,
                    entity_id: count
                };
                count += 1;

                let mut entity = PhyEntity::cube(
                    cube_id.clone(), Vector3::repeat(1.0)
                );
                entity.is.state.pos = Vector3::new(
                    x  as f64 * spacing - 5.0,
                    5.0 + y as f64 * spacing,
                    z as f64 * spacing - 5.0
                );
                entity.is.momentum = Vector3::new(0.0, -1.0, 0.0);
                entity.sync();

                commands
                    .spawn_bundle(PbrBundle {
                        mesh: cube_handle.clone(),
                        material: cube_material_handle.clone(),
                        transform: Transform::from_xyz(0.0, 0.0, 0.0),
                        ..default()
                    })
                    .insert(cube_id);

                engine.world.blas_mut().push(entity);
            }
        }
    }
    engine.world.build();


    unsafe {
        PhysicsEngine::init_global(engine);
    }



    // parent cube
    commands
        .spawn_bundle(PbrBundle {
            mesh: floor_handle.clone(),
            material: floor_material_handle.clone(),
            transform: Transform::from_xyz(0.0, 0.0, 0.0),
            .. default()
        })
        .insert(floor_id);

    // light
    commands.spawn_bundle(PointLightBundle {
        transform: Transform::from_xyz(4.0, 5.0, -4.0),
        .. default()
    });
    // camera
    commands.spawn_bundle(PerspectiveCameraBundle {
        transform: Transform::from_xyz(10.0, 15.0, 15.0)
            .looking_at(Vec3::new(0.0, 5.0, 0.0), Vec3::Y),
        .. default()
    });
}