(function () {
    'use strict';

    //const WIDTH = 1000;
    //const HEIGHT = 800;

    var projection = d3["geoOrthographic"]()
        .precision(0.05)
        .rotate([-92.9698998020521, -34.555785237703425, 0])
        .scale(500)
        .translate([500, 500])
        ;

    let sphere = { type: "Sphere" };
    let graticule = d3.geoGraticule10()

    var canvas = d3.select("canvas"),
        width = canvas.property("width"),
        height = canvas.property("height"),
        context = canvas.node().getContext("2d");

    const WIDTH = width;
    const HEIGHT = height;

    var timeSlider = document.getElementById('time_slider')

    var path = d3.geoPath()
        .projection(projection)
        .context(context);

    function renderBlank() {
        //context.clearRect(0, 0, WIDTH, HEIGHT);
        context.clearRect(0, 0, width, height);
        context.globalAlpha = 0.5;
        context.beginPath(), path(sphere), context.fillStyle = "#fff", context.fill();
        context.beginPath(), path(graticule), context.strokeStyle = "#ccc", context.stroke();
        context.beginPath(), path(sphere), context.stroke();

    }

    renderBlank()


    function attachEventListeners() {

        let drawBlocksButton = document.getElementById('draw-blocks')
        drawBlocksButton.addEventListener('click', drawBlocks)

    }



    function drawBlocks() {
        //let blocks_path = "./chn_blocks_simp.geojson"
        //let poles_path = "./block_poles_eur_rel.csv"

        //let range = (start, stop, step = 1) => Array(stop -
        //start).fill(start).map((x, y) => x + y * step)
        function range(start, stop, step = 1) {
            return Array(Math.ceil((stop - start) / step + 1)).fill(start).map((x, y) => x + y * step);
        }
        const blockTimes = range(parseFloat(timeSlider.min), parseFloat(timeSlider.max), parseFloat(timeSlider.step))

        Promise.all([
            d3.json(blocks_path),
            d3.csv(poles_path)
        ]).then(
            function (data) {
                // data[0].features.forEach(function (d, i) {
                //     d.geometry.coordinates.forEach(function (ring) {
                //         ring.reverse();
                //     });
                // });

                var poles = data[1].reduce(function (obj, x) {
                    obj[x.mov] = x;
                    return obj;
                }, {});

                var blocksAtTime = rotate_all_blocks(data[0], poles, blockTimes)

                //drawChart(geodata, poles);
                drawChart(blocksAtTime);

                timeSlider.addEventListener("input", function () {
                    drawChart(blocksAtTime)
                });
            });
    };



    function rotate_all_blocks(geodata_zero, poles, times) {

        let blocksAtTime = times.reduce(function (obj, time) {
            obj[time.toString()] = rotate_blocks(geodata_zero, poles, time);
            return obj;
        }, {});

        return blocksAtTime
    }

    function drawChart(blocksAtTime) {

        function render(blocksAtTime) {
            context.clearRect(0, 0, WIDTH, HEIGHT);
            context.globalAlpha = 0.5;
            context.beginPath(), path(sphere), context.fillStyle = "#fff", context.fill();
            context.beginPath(), path(graticule), context.strokeStyle = "#ccc", context.stroke();
            context.beginPath(), path(sphere), context.stroke();

            var geodata = blocksAtTime[timeSlider.value]

            geodata.features.forEach(function (d, i) {
                context.beginPath(), path(d), context.fillStyle = d.properties.color, context.fill();
                context.beginPath(), path(d), context.strokeStyle = d.properties.color, context.stroke();
            });
        }


        return d3.select(context.canvas)

            .call(zoom(projection)
                .on("zoom.render", () => render(blocksAtTime))
            )
            .call(() => render(blocksAtTime))
            .node();
    }


    function zoom(projection) {
        let v0, r0, q0;

        function zoomstarted() {
            v0 = versor.cartesian(projection.invert(d3.mouse(this)));
            r0 = projection.rotate();
            q0 = versor(r0);
        }

        function zoomed() {
            projection.scale(d3.event.transform.k * (height - 10) / 2);

            var v1 = versor.cartesian(projection.rotate(r0).invert(d3.mouse(this))),
                q1 = versor.multiply(q0, versor.delta(v0, v1)),
                r1 = versor.rotation(q1);
            projection.rotate(r1);
        }

        return d3.zoom()
            .on("start", zoomstarted)
            .on("zoom", zoomed)
    }


    attachEventListeners();
})();
