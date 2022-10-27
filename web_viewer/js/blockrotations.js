function radian(deg) {
    return deg * Math.PI / 180
}

function degree(rad) {
    return rad * 180 / Math.PI
}

function sind(ang) {
    return Math.sin(radian(ang))
}

function cosd(ang) {
    return Math.cos(radian(ang))
}

function atand(x) {
    return degree(Math.atan(x))
}

function atand2(y, x) {
    return degree(Math.atan2(y, x))
}

function sphere_to_cart(lon, lat) {

    x = cosd(lat) * cosd(lon)
    y = cosd(lat) * sind(lon)
    z = sind(lat)

    return [x, y, z]
}

function cart_to_sphere(x, y, z) {

    let lon = atand2(y, x)
    let lat = atand(z / Math.sqrt(x * x + y * y))

    return [lon, lat]
}


function coord_sphere_to_cart(lonlat) {
    return sphere_to_cart(lonlat[0], lonlat[1])
}


function coord_cart_to_sphere(xyz) {
    return cart_to_sphere(xyz[0], xyz[1], xyz[2])
}


function make_rot_matrix(cpole, ang) {

    let Ex = cpole[0]
    let Ey = cpole[1]
    let Ez = cpole[2]

    let r11 = Ex * Ex * (1 - cosd(ang)) + cosd(ang)
    let r12 = Ex * Ey * (1 - cosd(ang)) - Ez * sind(ang)
    let r13 = Ex * Ez * (1 - cosd(ang)) + Ey * sind(ang)

    let r21 = Ey * Ex * (1 - cosd(ang)) + Ez * sind(ang)
    let r22 = Ey * Ey * (1 - cosd(ang)) + cosd(ang)
    let r23 = Ey * Ez * (1 - cosd(ang)) - Ex * sind(ang)

    let r31 = Ez * Ex * (1 - cosd(ang)) - Ey * sind(ang)
    let r32 = Ez * Ey * (1 - cosd(ang)) + Ex * sind(ang)
    let r33 = Ez * Ez * (1 - cosd(ang)) + cosd(ang)

    let row1 = [r11, r12, r13]
    let row2 = [r21, r22, r23]
    let row3 = [r31, r32, r33]

    let rotMat = [row1, row2, row3]

    return rotMat
}


function matrix_dot(A, B) {
    let result = new Array(A.length).fill(0).map(row => new Array(B[0].length).fill(0));

    return result.map((row, i) => {
        return row.map((val, j) => {
            return A[i].reduce((sum, elm, k) => sum + (elm * B[k][j]), 0)
        })
    })
}

function multiply(a, b) {
    var aNumRows = a.length, aNumCols = a[0].length,
        bNumRows = b.length, bNumCols = b[0].length,
        m = new Array(aNumRows);  // initialize array of rows
    for (var r = 0; r < aNumRows; ++r) {
        m[r] = new Array(bNumCols); // initialize the current row
        for (var c = 0; c < bNumCols; ++c) {
            m[r][c] = 0;             // initialize the current cell
            for (var i = 0; i < aNumCols; ++i) {
                m[r][c] += a[r][i] * b[i][c];
            }
        }
    }

    return m;
}

function matrixProduct(matrix1, matrix2) {
    // Assumes matrix1 and matrix 2 can be multiplied
    // but one can check if the dimensions match here 

    // We know that for matrices [m x k] X [k x n] 
    // that the result is [m x n] and the k's are 
    // the dot product matching lengths
    let result = [] // if for loops are correct then the dims will be correct 

    // For intuition, m is the number of rows of matrix1
    // n is the number of columns of matrix2
    // And in our 2D arrays, number of rows is matrix.length (outer array)
    // number of columns is matrix[i].length (inner array or nested array)

    // For consistency, I will call X the index of row of matrix 1
    // and Y the index of column of matrix 2 
    for (let X = 0; X < matrix1.length; X++) {

        let newRow = [] // new row to add to result, and to add scalars to

        for (let Y = 0; Y < matrix2[0].length; Y++) {
            // Once we've "arrived" at position X, Y
            // we want to find the dot product of row X of matrix 1
            // and column Y of matrix 2 (which are the same length) 

            let newDotProd = 0 // add together products to get dot product

            for (let i = 0; i < matrix1[X].length; i++) {
                // note the switched positions of the i index
                newDotProd += matrix1[X][i] * matrix2[i][Y]

            }

            newRow.push(newDotProd) // should have pushed n times
        }

        result.push(newRow) // should have pushed m times
    }

    return result
} // end of matrixProduct()



function copy(o) {
    var output, v, key;
    output = Array.isArray(o) ? [] : {};
    for (key in o) {
        v = o[key];
        output[key] = (typeof v === "object") ? copy(v) : v;
    }
    return output;
}


function rotate_block_coords(coords, pole, time) {

    var cpole = sphere_to_cart(pole["lon"], pole["lat"])
    var rot_matrix = make_rot_matrix(cpole, pole["rotrate"] * time)

    var cart_trace = [coords[0].map(coord_sphere_to_cart)]

    var new_cart_coords = [cart_trace[0].map(function (coord) {
        return matrix_dot([coord], rot_matrix)[0]
    })]

    var new_coords = [new_cart_coords[0].map(coord_cart_to_sphere)]

    // debugger;

    return new_coords
}

function rotate_block(block, pole, time) {
    block.geometry.coordinates = rotate_block_coords(
        block.geometry.coordinates, pole, time
    )
}

function rotate_blocks(blocks, poles, time) {
    var new_blocks = copy(blocks)

    new_blocks.features.forEach(function (d) {
        var pole = poles[d.properties.fid.toString()];
        rotate_block(d, pole, time)
    }
    )

    return new_blocks
}