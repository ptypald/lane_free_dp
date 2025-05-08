// var scene_el = document.getElementById("myDiv");
// var road = scene_el.getBoundingClientRect();
// var roadw = road.width;
// // var roadw = window.innerWidth;
// var roadh = road.height

var conf = {
    el: "#app",
    data: {
        sim: window.sim,
        // t: 0,
        px_per_meter: 1,
		px_per_meter_max: 10,
		px_per_meter_min: 1,
		d3: {
			data: null,
            data_crashes: null,
            data_progress: 0,
		},
        // d3_data_progress: 0,
        play: {
			t: 0,
			active: false,
			intvl: null,
			speed: 3
        },
        chart: {
			mode: "general",
			vehicle: {
				id: 0
			},
			general: {
				type: "mdsd"
			}
		},
        vehSelect: {
            active: false,
            id: null,
            text: '',
            marker: false
        },
        // vehicle Theme
        rafReference: null,   //check
        /* toolbar */
        switchTheme: true,
        showDropMenu: true,
        items:[
            { title: 'FAQs' },
            { title: 'About' }
        ],
        logos:[
            {title: 'DSSL', url: 'http://www.dssl.tuc.gr'},
            {title: 'TUC', url:"http://www.tuc.gr/"},
            {title: 'TRAFFICFLUID', url:"https://www.trafficfluid.tuc.gr/"}
        ],

        myThreeView: {
            pov: false,
            top: false
        }

    },
    computed: {
        trackTime(){
            if (this.vehSelect.id != null)
                this.vehSelect.active = true;
            // dummy time step variable for watch: {}
            return this.play.t;
        },
    },

    watch: {
        // whenever trackTime changes, this function will run
        trackTime() {
            // console.log(`${this.play.t}`);
            // this.svg_show(this.play.t);
        },

    },

    methods: {
        
        road_container_resize(el) {
            let width_px = $(el).width() - 10;
            this.px_per_meter = Math.min(this.px_per_meter_max, Math.max(this.px_per_meter_min, width_px / this.sim.roadlength));
            if (!this.d3.svg_main)
                this.svg_initialize();
            // else
            // this.threeEnvironment();
            // this.svg_show(this.play.t);
        },

        svg_initialize() {
            var scene_el = document.getElementById("myDiv");
            var road = scene_el.getBoundingClientRect();
            var roadw = road.width;
            var roadh = road.height

            this.d3.svg_main = d3.select("#myDiv")
                .attr("width", roadw)
                .attr("height", roadh);
            this.svg_load_data();
            this.threeEnvironment();

        },

        threeEnvironment() {
            var scene_el = document.getElementById("myDiv");
            var road = scene_el.getBoundingClientRect();
            var roadw = road.width;
            var roadh = road.height
            
            var pixelRatio = window.devicePixelRatio;
            var sceneWidth = roadw * pixelRatio;
            var sceneHeight = roadh * pixelRatio;

            const scene = new THREE.Scene();
            scene.background = new THREE.Color( 0xcccccc );

            raycaster = new THREE.Raycaster();

            const renderer = new THREE.WebGLRenderer({
                scene_el,
                alpha: true,
                antialias: true,
                // powerPreference: "high-performance"
            });
            renderer.setPixelRatio( pixelRatio );
            renderer.shadowMap.enabled = true;
            renderer.shadowMap.type = THREE.PCFSoftShadowMap;
            // renderer.setSize( window.innerWidth, window.innerHeight );
            renderer.setSize( sceneWidth, sceneHeight );
            // document.body.appendChild( renderer.domElement );
            document.getElementById("myDiv").appendChild( renderer.domElement );

            var allVeh = [];
            const frustumSize = 25;
            // const viewWidth = window.innerWidth
            // const viewHeight = window.innerHeight
            // const aspect = window.innerWidth / window.innerHeight;
            const viewWidth = sceneWidth
            const viewHeight = sceneHeight
            const aspect = sceneWidth / sceneHeight;
            const camera = new THREE.OrthographicCamera( frustumSize * aspect / - 2, frustumSize * aspect / 2, frustumSize / 2, frustumSize / - 2, 0.1, 10000 );

            // const camera = new THREE.PerspectiveCamera( 160, window.innerWidth / window.innerHeight, 1, 1000 );
            camera.position.set( 0, -5, 10 );

            var zoom = 0.1;

            /* ADD for camera controls zoom in/out */
            const controls = new THREE.OrbitControls( camera, renderer.domElement );
            controls.enableRotate = false;
            // controls.listenToKeyEvents( window ); // optional
            controls.enableDamping = true; // an animation loop is required when either damping or auto-rotation are enabled
            controls.dampingFactor = 0.05;
            controls.screenSpacePanning = true;
            controls.maxZoom = 10.0;
            controls.minZoom = 0.4;

            controls.mouseButtons = {
                LEFT: THREE.MOUSE.PAN,
                MIDDLE: THREE.MOUSE.DOLLY,
                RIGHT: THREE.MOUSE.ROTATE
            }

            let lanes;

            const carFrontTexture = new Texture(40,80,[{x: 0, y: 10, w: 30, h: 60 }]);
            const carBackTexture = new Texture(40,80,[{x: 10, y: 10, w: 30, h: 60 }]);
            const carRightSideTexture = new Texture(110,40,[{x: 10, y: 0, w: 30, h: 30 }, {x: 50, y: 0, w: 50, h: 30 }]);
            const carLeftSideTexture = new Texture(110,40,[{x: 10, y: 10, w: 30, h: 30 }, {x: 50, y: 10, w: 50, h: 30 }]);

            // const generateLanes = () => [0, 1, 2].map((index) => {
            const generateLanes = () => [0, 1].map((index) => {
                const lane = new Lane(index);
                
                // lane.mesh.position.y = index*positionWidth*zoom;
                scene.add( lane.mesh );

                scene.add( lane.mesh );

                return lane;
            }).filter((lane) => lane.index >= 0);

            hemiLight = new THREE.HemisphereLight(0xffffff, 0xffffff, 0.6);
            scene.add(hemiLight)

            const initialDirLightPositionX = -100;
            const initialDirLightPositionY = -100;
            dirLight = new THREE.DirectionalLight(0xffffff, 0.6);
            dirLight.position.set(initialDirLightPositionX, initialDirLightPositionY, 200);
            dirLight.castShadow = true;
            // dirLight.target = chicken;
            scene.add(dirLight);

            // dirLight.shadow.mapSize.width = window.innerWidth;
            // dirLight.shadow.mapSize.height = window.innerHeight;
            dirLight.shadow.mapSize.width = sceneWidth;
            dirLight.shadow.mapSize.height = sceneHeight;
            var d = 100;
            dirLight.shadow.camera.left = - d;
            dirLight.shadow.camera.right = d;
            dirLight.shadow.camera.top = d;
            dirLight.shadow.camera.bottom = - d;

            // var helper = new THREE.CameraHelper( dirLight.shadow.camera );
            // var helper = new THREE.CameraHelper( camera );
            // scene.add(helper)

            backLight = new THREE.DirectionalLight(0x000000, .4);
            backLight.position.set(200, 200, 50);
            backLight.castShadow = true;
            scene.add(backLight)

            const laneTypes = ['forest', 'car', 'boundary'];
            const vechicleColors = [0xa52523, 0xbdb638, 0x78b14b];
            const threeHeights = [2.0, 3.5,4.0];

            const initaliseValues = () => {
                lanes = generateLanes()

                dirLight.position.x = initialDirLightPositionX;
                dirLight.position.y = initialDirLightPositionY;

                window.addEventListener( 'pointermove', onPointerMove );
            }

            initaliseValues();

            function Texture(width, height, rects) {
                const canvas = document.createElement( "canvas" );
                canvas.width = width;
                canvas.height = height;
                const context = canvas.getContext( "2d" );
                context.fillStyle = "#ffffff";
                context.fillRect( 0, 0, width, height );
                context.fillStyle = "rgba(0,0,0,0.6)";  
                rects.forEach(rect => {
                    context.fillRect(rect.x, rect.y, rect.w, rect.h);
                });

                return new THREE.CanvasTexture(canvas);
            }

            function Wheel(id) {
                const wheel = new THREE.Mesh( 
                    new THREE.BoxBufferGeometry( (sim.Cx[id]/4)*zoom, (sim.Cy[id] / 3)*zoom, 1.2*zoom ), 
                    new THREE.MeshPhongMaterial( { color: 0x333333, flatShading: true } ) 
                );
                wheel.position.z = 0.6*zoom;
                return wheel;
            }

            function Car(id) {
                const car = new THREE.Group();
                var color;

                if (sim.id[id] == 1) {
                    color = 0xD80621;
                }
                else if (sim.id[id] == 2) {
                    color = 0x1493ff;
                    // color = 0x4C00FF;
                }
                else
                    color = vechicleColors[1];    
                    // color = vechicleColors[Math.floor(Math.random() * vechicleColors.length)];

                const main = new THREE.Mesh(
                    new THREE.BoxBufferGeometry( sim.Cx[id]*zoom, sim.Cy[id]*zoom, 1.5*zoom ), 
                    new THREE.MeshPhongMaterial( { color, flatShading: true } )
                );
                main.position.z = 1.2*zoom;
                main.castShadow = true;
                main.receiveShadow = true;
                car.add(main)

                const cabin = new THREE.Mesh(
                    new THREE.BoxBufferGeometry( sim.Cx[id]*zoom*0.7, sim.Cy[id]*zoom*0.8, 1.2*zoom ), 
                    [
                        new THREE.MeshPhongMaterial( { color: 0xcccccc, flatShading: true, map: carBackTexture } ),
                        new THREE.MeshPhongMaterial( { color: 0xcccccc, flatShading: true, map: carFrontTexture } ),
                        new THREE.MeshPhongMaterial( { color: 0xcccccc, flatShading: true, map: carRightSideTexture } ),
                        new THREE.MeshPhongMaterial( { color: 0xcccccc, flatShading: true, map: carLeftSideTexture } ),
                        new THREE.MeshPhongMaterial( { color: 0xcccccc, flatShading: true } ), // top
                        new THREE.MeshPhongMaterial( { color: 0xcccccc, flatShading: true } ) // bottom
                    ]
                );

                
                cabin.position.x = -0.2*zoom;
                cabin.position.z = 2.55*zoom;
                cabin.castShadow = true;
                cabin.receiveShadow = true;
                car.add( cabin );

                const frontRightWheel = new Wheel(id);
                frontRightWheel.position.x = -(sim.Cx[id]/3)*zoom;
                frontRightWheel.position.y = -(sim.Cy[id]/3 + 0.1)*zoom;
                car.add( frontRightWheel );

                const frontLeftWheel = new Wheel(id);
                frontLeftWheel.position.x = -(sim.Cx[id]/3)*zoom;
                frontLeftWheel.position.y = (sim.Cy[id]/3 + 0.1)*zoom;
                car.add( frontLeftWheel );

                const backRightWheel = new Wheel(id);
                backRightWheel.position.x = (sim.Cx[id]/3)*zoom;
                backRightWheel.position.y = -(sim.Cy[id]/3 + 0.1)*zoom;
                car.add( backRightWheel );

                const backLeftWheel = new Wheel(id);
                backLeftWheel.position.x = (sim.Cx[id]/3)*zoom;
                backLeftWheel.position.y = (sim.Cy[id]/3 + 0.1)*zoom;
                car.add( backLeftWheel );

                car.castShadow = true;
                car.receiveShadow = false;

                return car;  
            }

            function Road() {
                const road = new THREE.Group();

                const createSection = color => new THREE.Mesh(
                    // new THREE.PlaneBufferGeometry( sim.roadlength*(zoom), sim.roadwidth*(zoom) ),
                    new THREE.BoxBufferGeometry( sim.roadlength*zoom, sim.roadwidth*zoom, 0.01*zoom ), 
                    new THREE.MeshPhongMaterial( { color } )
                );

                const middle = createSection(0x454A59);
                middle.receiveShadow = true;
                // middle.position.y = 100*zoom
                middle.position.x = (sim.roadlength/2.0)*zoom;
                road.add(middle);

                return road;
            }

            function SelectedMarker() {

                const marker = new THREE.Group();

                const createSection = color => new THREE.Mesh(
                    // new THREE.PlaneBufferGeometry( 5.0*(zoom), 2.0*(zoom) ),
                    new THREE.BoxBufferGeometry( 10.0*zoom, 2.0*zoom, 0.1*zoom ), 
                    new THREE.MeshPhongMaterial( { color } )
                );

                const groundMarker = createSection(0xbaf455);
                groundMarker.receiveShadow = true;
                // middle.position.y = 100*zoom
                // middle.position.x = (sim.roadlength/2.0)*zoom;
                marker.add(groundMarker);

                return marker;
            }

            function Ramp() {
                
            }

            function Boundary () {
                var sections = sim.offset_count + 1;
                var lim1, lim2, slope, offset, x;        
                var y = {};
                var length = 1500;

                for (var i = 0; i < length; i++) {
                    var temp_y = 0
                    for (j = 0; j < sections - 1; j++) {
                        lim1   = sim.limits[j][i];
                        lim2   = sim.limits[j+1][i];
                        slope  = sim.slopes[j];
                        offset = sim.offsets[j];
                        x = i;
                        temp_y += 0.5 * (lim2 - lim1) * Math.tanh(slope * (x - offset));
                    }
                    y[i] = sim.roadwidth - (temp_y + 0.5*(sim.limits[0][0] + sim.limits[sections - 1][0]));
                }
                const blocks = new THREE.Group();

                const createSection = (color, width) => new THREE.Mesh(
                    new THREE.BoxBufferGeometry( 1*zoom, (0.1 + width)*zoom, 1.0*zoom ), 
                    new THREE.MeshPhongMaterial( { color } )
                );

                var block = [];
                for (var i = 0; i < length; i++) {
                    block[i] = createSection(0xbaf455, y[i]);
                    block[i].receiveShadow = true;
                    block[i].position.x = ((i + 1) *1.5 )*zoom;
                    block[i].position.y = ((sim.roadwidth - y[i]) / 2.0) *zoom;
                    blocks.add(block[i]);
                }

                return blocks;

            }

            function Grass() {
                const grass = new THREE.Group();

                const createSection = color => new THREE.Mesh(
                    new THREE.BoxBufferGeometry( sim.roadlength*zoom, sim.roadwidth*zoom, 1.0*zoom ), 
                    new THREE.MeshPhongMaterial( { color } )
                );

                const right = createSection(0xbaf455);
                right.receiveShadow = true;
                right.position.x = (sim.roadlength/2.0)*zoom;
                right.position.y = - sim.roadwidth*zoom;
                grass.add(right);

                const left = createSection(0xbaf455);
                left.receiveShadow = true;
                left.position.x = (sim.roadlength/2.0)*zoom;
                left.position.y = sim.roadwidth*zoom;
                grass.add(left);

                // grass.position.z = 1.5*zoom;
                return grass;
            }

            function Tree() {
                const three = new THREE.Group();
            
                const trunk = new THREE.Mesh(
                new THREE.BoxBufferGeometry( 1.5*zoom, 1.5*zoom, 4*zoom ), 
                new THREE.MeshPhongMaterial( { color: 0x4d2926, flatShading: true } )
                );
                trunk.position.z = 2*zoom;
                trunk.castShadow = true;
                trunk.receiveShadow = true;
                three.add(trunk);
            
                height = threeHeights[Math.floor(Math.random()*threeHeights.length)];
            
                const crown = new THREE.Mesh(
                new THREE.BoxBufferGeometry( 3*zoom, 3*zoom, height*zoom ), 
                new THREE.MeshPhongMaterial( { color: 0x7aa21d, flatShading: true } )
                );
                crown.position.z = (height/2+4)*zoom;
                crown.castShadow = true;
                crown.receiveShadow = false;
                three.add(crown);
            
                return three;  
            }

            function Lane(index) {
                this.index = index;
                this.type = index <= 0 ? 'forest' : laneTypes[index];

                switch(this.type) {
                    case 'field': {
                        this.type = 'field';
                        this.mesh = new Grass();
                        break;
                    }
                    case 'boundary': {
                        this.type = 'boundary';
                        this.mesh = new Boundary();
                        break;
                    }
                    case 'forest': {
                        this.mesh = new Grass();
                        
                        var treeSpacing = 200.0;

                        var numTrees = Math.floor(sim.roadlength / treeSpacing) - 1;
                        numTrees *= 2;
                        treeList = Array(numTrees).fill().map((_, idx) => 1 + idx);
                        
                        this.occupiedPositions = new Set();
                        this.trees = treeList.map((index) => {
                            const tree = new Tree();
                        
                            treePos = index % (numTrees/2) + 1;
                            tree.position.x = treePos * treeSpacing * zoom;
                            
                            tree.position.y = (sim.roadwidth) * zoom;
                            if (index > numTrees / 2) {
                                tree.position.y = -(sim.roadwidth) * zoom;
                            }
                            this.mesh.add( tree );
                            
                            return tree;
                        })
                        break;
                    }
                    case 'car' : {
                        this.mesh = new Road();

                        let A = [...Array(sim.n).keys()];;
                        this.vechicles = A.map((el, id) => {

                            // here add ambulance and police vehicle
                            const vechicle = new Car(id);

                            vechicle.position.x = sim.x[id][0]*zoom;
                            vechicle.position.y = (sim.y[id][0] - (sim.roadwidth / 2.0)) * zoom;
                            // vechicle.position.y = (sim.y[id][0]) * zoom;
                            
                            this.mesh.add( vechicle );
                            allVeh.push( vechicle );
                            
                            return vechicle;
                        })

                        break;
                    }
                }
            }

            var vec3 = new THREE.Vector3();
            const pointer = new THREE.Vector2();
            let INTERSECTED;

            function onPointerMove( event ) {
                // calculate pointer position in normalized device coordinates
                // (-1 to +1) for both components
                pointer.x = ( event.clientX / sceneWidth ) * 2 - 1;
                pointer.y = - ( event.clientY / sceneHeight ) * 2 + 1;
            }
            
            function animate() {    
                requestAnimationFrame( animate );
                // console.log("run");
                count_time = conf.data.play.t;
                var tweenTime = 0.125; // conf.data.sim.Step; // (conf.data.sim.Step*1000/conf.data.play.speed);
                if ( count_time < sim.k) {
                    // Animate cars and trucks moving on the lane
                    lanes.forEach(lane => {
                        if(lane.type === 'car') {
                            lane.vechicles.forEach((vechicle, id )=> {
                                vec3.subVectors(camera.position, vechicle.position);
                                
                                if (count_time == 0) {
                                    vechicle.position.x = (sim.x[id][count_time]*zoom);
                                    vechicle.position.y = (sim.y[id][count_time] - (sim.roadwidth / 2.0)) * zoom ;
                                }
                                else if ((sim.x[id][count_time] < sim.x[id][count_time - 1])) {
                                    // vechicle.position.x = (sim.x[id][count_time]*zoom);
                                    // vechicle.position.y = (sim.y[id][count_time] - (sim.roadwidth / 2.0)) * zoom ;
                                    var tween = new TWEEN.Tween(vechicle.position) // Create a new tween that modifies 'coords'.
                                        .to({x: sim.x[id][count_time]*zoom, y: (sim.y[id][count_time] - (sim.roadwidth / 2.0)) * zoom}, 0.0) // Move to (300, 200) in 1 second.
                                        .easing(TWEEN.Easing.Linear.None) // Use an easing function to make the animation smooth.
                                        .onUpdate(() => {
                                            // console.log("hello");
                                            // Called after tween.js updates 'coords'.
                                            // Move 'box' to the position described by 'coords' with a CSS translation.
                                            // box.style.setProperty('transform', `translate(${coords.x}px, ${coords.y}px)`)
                                        })
                                        .start() // Start the tween immediately.
                                }
                                else {
                                    var tween = new TWEEN.Tween(vechicle.position) // Create a new tween that modifies 'coords'.
                                        .to({x: sim.x[id][count_time]*zoom, y: (sim.y[id][count_time] - (sim.roadwidth / 2.0)) * zoom}, tweenTime) // Move to (300, 200) in 1 second.
                                        .easing(TWEEN.Easing.Linear.None) // Use an easing function to make the animation smooth.
                                        .onUpdate(() => {
                                            // console.log("hello");
                                            // Called after tween.js updates 'coords'.
                                            // Move 'box' to the position described by 'coords' with a CSS translation.
                                            // box.style.setProperty('transform', `translate(${coords.x}px, ${coords.y}px)`)
                                        })
                                        .start() // Start the tween immediately.
                                }

                                // if follows veh with specific ID
                                // if (conf.data.vehSelect.id != null && conf.data.vehSelect.id == id) {
                                if (conf.data.vehSelect.id != null && conf.data.vehSelect.id == sim.id[id]) {
                                    
                                    controls.object.position.copy(vechicle.position).add(vec3);
                                    controls.target.copy(vechicle.position);
        
                                    if (!conf.data.myThreeView.pov && !conf.data.myThreeView.top) {
                                        camera.position.x = vechicle.position.x;
                                        camera.position.y = -5;
                                        camera.position.z = 10;
                                    }
                                    else if (conf.data.myThreeView.top) {
                                        camera.position.x = vechicle.position.x;
                                        camera.position.y = 0;
                                        camera.position.z = 50;
                                    }
                                    else if (conf.data.myThreeView.pov) {
                                        // TODO
                                    }

                                    if (conf.data.vehSelect.marker == false) {
                                        conf.data.vehSelect.marker = true
                                        var myMarker = new SelectedMarker();
                                        
                                        myMarker.position = vechicle.position;                                        
                                        myMarker.name = "myMarker";
                                        scene.add( myMarker );
                                    }
                                    
                                }
                                else if (conf.data.vehSelect.id == null && conf.data.vehSelect.marker == true){
                                    // console.log("remove marker");
                                    var test = scene.getObjectByName("myMarker");
                                    // console.log(test);
                                    test.parent.remove(test);
                                    // test.geometry.dispose();

                                    conf.data.vehSelect.marker = false
                                }
                            });
                        }
                    });
                }
                controls.update();
                // TWEEN.update();
                renderer.render( scene, camera );
                TWEEN.update();
            }

            const resizeUpdateInterval = 500;

            window.addEventListener(
            'resize',
            // throttle(
                () => {
                    var scene_el = document.getElementById("myDiv");
                    var road = scene_el.getBoundingClientRect();
                    var roadw = road.width;
                    var roadh = road.height

                    camera.aspect = roadw / roadh;
                    camera.updateProjectionMatrix();
                    renderer.setSize(roadw, roadh);
                },
            //     resizeUpdateInterval,
            //     { trailing: true }
            // )
            );

            animate();
        },

        svg_load_data() {
            // create d3 data vector
            this.d3.data = [];
            for (var t = 0; t < this.sim.k; t++) {
                this.d3.data_progress = Math.round(((t/this.sim.k) * 100));
                this.d3.data.push([]);
                var j = 0;
                for (var i=0; i < this.sim.n; i++) {
                    // create data[t][i] = {}
                    this.d3.data[t].push({});
                    // this.d3.data[t][j].id = j;
                    this.d3.data[t][j].id = this.sim.id[j];
                    this.d3.data[t][j].l = this.sim.Cx[j];
                    this.d3.data[t][j].w = this.sim.Cy[j];
                    // this.d3.data[t][j].vdx = this.sim.vdx[j];                    
                    this.d3.data[t][j].x = this.sim.x[j][t];
                    this.d3.data[t][j].y = this.sim.roadwidth - this.sim.y[j][t];
                    // this.d3.data[t][j].vx = this.sim.vx[j][t];

                    // calculate veh. angle data
                    var temp_angle;                    
                    if (t != 0) {
                        if (this.sim.y[j][t - 1] == this.sim.y[j][t])
                            temp_angle = 0;
                        else
                            temp_angle = Math.atan2(this.sim.y[j][t] - this.sim.y[j][t - 1], this.sim.x[j][t] - this.sim.x[j][t - 1]);
                    }
                    else
                        temp_angle = 0;                                        
                    this.d3.data[t][j].angle = -temp_angle * (180/Math.PI);
                    
                    j++;
                }
            }
            this.d3.data_progress = 100;
        },
        
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//        
        play_start() {
            this.play.active = true;
            // this.svg_show(this.play.t);
            this.threeEnvironment;
            
            var fpsTime = (this.sim.Step*1000/this.play.speed);

            this.play.intvl = setInterval(() => {
                if (this.play.t >= this.sim.k-1) {
                    this.play_stop();
                    this.play.t = 0;
                    this.vehSelect.active = false;
                    this.vehSelect.id = null;
                    return;
                }
                this.play.t += 1;
                // this.svg_show(this.play.t);
                this.threeEnvironment;
                
            }, fpsTime);
        },

        play_back_start() {
            this.play.active = true;
            // this.svg_show(this.play.t);
            this.threeEnvironment;
            
            var fpsTime = (this.sim.Step*1000/this.play.speed);

            this.play.intvl = setInterval(() => {
                if (this.play.t >= this.sim.k-1) {
                    this.play_stop();
                    this.play.t = 0;
                    this.vehSelect.active = false;
                    this.vehSelect.id = null;
                    return;
                }
                this.play.t -= 1;
                // this.svg_show(this.play.t);
                this.threeEnvironment;
                
            }, fpsTime);
        },

        play_stop() {
            this.play.active = false;            
            clearInterval(this.play.intvl);
        },

        play_next_step() {
            return this.play.t += 1;
        },

        play_prev_step() {
            return this.play.t -= 1;
        },

        chart_show() {
            if (this.vehSelect.active)
                window.chart_show(this);
            else
                console.log("No vehicle selected");
        }
        
    },
};

new Vue(window.conf);
Vue.use(Vuetify);

