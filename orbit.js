
/****************************************
 * Kepler's equation solvers 
 * Mostly copied, with modifications,
 * from http://www.jgiesen.de/kepler/kepler.html
 ****************************************/

 function EccAnom(ec,m,dp) {

    // arguments:
    // ec=eccentricity, m=mean anomaly,
    // dp=number of decimal places
    var pi=Math.PI, K=pi/180.0;

    var maxIter=30, i=0;

    var delta=Math.pow(10,-dp);

    var E, F;

    m=m/360.0;

    m=2.0*pi*(m-Math.floor(m));

    if (ec<0.8) E=m; else E=pi;

    F = E - ec*Math.sin(m) - m;

    while ((Math.abs(F)>delta) && (i<maxIter)) {

      E = E - F/(1.0-ec*Math.cos(E));
      F = E - ec*Math.sin(E) - m;

      i = i + 1;

    }

    E=E/K;

    return E;

  }

  function TrueAnom(ec,E) {

    K=Math.PI/180.0;
    S=Math.sin(E);

    C=Math.cos(E);

    fak=Math.sqrt(1.0-ec*ec);

    phi=Math.atan2(fak*S,C-ec)/K;

    return phi;

  }


  function position(a, ec, E) {

    // a=semimajor axis, ec=eccentricity, E=eccentric anomaly
    // x,y = coordinates of the planet with respect to the Sun

    C = Math.cos(E);


    S = Math.sin(E);

    x = a*(C-ec);

    y = a*Math.sqrt(1.0-ec*ec)*S;

    return [x, y];

  }


/*****************************
* Math utilities
******************************/ 
function sign (x) {
  return x != 0 ? Math.abs(x) / x : 0; 
};

/*****************************
* This plots and does everything else
******************************/ 
doAll();
function doAll(){

    var RADFROMDEG = Math.PI / 180.;
    var TWOPI = 2. * Math.PI; 

    /**********************
    * Plotting options 
    ***********************/

    var width = 640;
    var height = 330;
    var plotHeight = height / 2;
    var padding = 30;

    var ratio = height / width;
    var plotRatio = plotHeight / width;

    /***********************
    * Physical quantities 
    ************************/

    var pulsePeriod = 0.1;
    var orbitPeriod = 6;

    var omega =  parseFloat(document.getElementById("omegaInput").value);
    var eccentricity = parseFloat(document.getElementById("eccInput").value);

    if (isNaN(omega) || isNaN(eccentricity) ||
       (eccentricity < 0) || (eccentricity >= 1)){
      console.log("Entered values are wrong")
      return -1;
    };

    var semimajAxis = 150; // In pixels, for now
    var tMax = 10000;
    var tMin = 0;
    var tSpan = 10000;

    var minPeriod = pulsePeriod;
    var maxPeriod = pulsePeriod;
    var minDelay = 0;
    var maxDelay = 0;

    // Derived quantities
    var omegaPSky = omega + 180; // Omega is calculated from the plane of the sky.

    if (eccentricity == 0) {omegaPSky = 180;};

    var omegaPSkyRAD = omegaPSky * RADFROMDEG;

    var semiminAxis = Math.sqrt(semimajAxis * semimajAxis * (1 - eccentricity * eccentricity));
    var focusDist = Math.sqrt(semimajAxis * semimajAxis - semiminAxis * semiminAxis);
    var focusOffX = focusDist * Math.cos(omegaPSkyRAD);
    var focusOffY = focusDist * Math.sin(omegaPSkyRAD);
    console.log(omega, omegaPSky, eccentricity, semimajAxis, semiminAxis);

    var phi0_x = semimajAxis * Math.cos(omegaPSkyRAD) + (width/2); 
    var phi0_y = semimajAxis * Math.sin(omegaPSkyRAD) + (height/2); 

    var phi90_x = semiminAxis * Math.cos(omegaPSkyRAD + Math.PI * 0.5) + (width/2); 
    var phi90_y = semiminAxis * Math.sin(omegaPSkyRAD + Math.PI * 0.5) + (height/2); 

    var phi180_x = semimajAxis * Math.cos(omegaPSkyRAD + Math.PI) + (width/2); 
    var phi180_y = semimajAxis * Math.sin(omegaPSkyRAD + Math.PI) + (height/2); 

    var phi270_x = semiminAxis * Math.cos(omegaPSkyRAD + Math.PI * 1.5) + (width/2); 
    var phi270_y = semiminAxis * Math.sin(omegaPSkyRAD + Math.PI * 1.5) + (height/2); 

    var distANode_x = semimajAxis * (1 - eccentricity * eccentricity) / (1 - eccentricity * Math.cos(omegaPSkyRAD))
    var distDNode_x = semimajAxis * (1 - eccentricity * eccentricity) / (1 + eccentricity * Math.cos(omegaPSkyRAD))

    /**********************
    * Define SVG canvas 
    ***********************/

    var svg = d3.select("#coolStuff").append("svg")
                              .attr("width", width)
                              .attr("height", height)
                              .append("g")
                              .attr("class", "orbitSVG");

    var svgYScale = d3.scale.linear()
                            .domain([-100, 100])
                            .range([height - padding, padding]); // even better with pa

    var svgXScale = d3.scale.linear()
                            .domain([-100 / ratio, 100 / ratio])
                            .range([padding, width - padding]); // even better with pa

    var svgPeriods = d3.select("#coolStuff").append("svg")
                                      .attr("width", width)
                                      .attr("height", plotHeight)
                                      .append("g")
                                      .attr("class", "orbitSVG");

    var svgPerYScale = d3.scale.linear()
                               .domain([minPeriod, maxPeriod])
                               .range([plotHeight - padding, padding]); 

    var svgPerXScale = d3.scale.linear()
                               .domain([0, tMax])
                               .range([padding, width - padding]); 

    var svgPerXAxis = d3.svg.axis()
                            .scale(svgPerXScale)
                            .orient("bottom") // orientation of labels
                            .ticks(5);  //Set rough # of ticks. D3 will take this as suggestion and decide based on "nice", whole intervals

    svgPeriods.append("g")
        .attr("class", "xaxis")  //Assign "axis" class
        .attr("transform", "translate(0," + (plotHeight - padding) + ")") // To align it to the bottom
        .call(svgPerXAxis);

    var svgPerYAxis = d3.svg.axis()
        .scale(svgPerYScale)
        .orient("left") // orientation of labels
        .ticks(5);  //Set rough # of ticks. D3 will take this as suggestion and decide based on "nice", whole intervals

    svgPeriods.append("g")
        .attr("class", "yaxis")
        .attr("transform", "translate(" + padding + ",0)")
        .call(svgPerYAxis);

    svgPeriods.append("text")
        .attr("class", "y label")
        .attr("text-anchor", "end")
        .attr("y", 6)
        .attr("x", -20)
        .attr("dy", ".75em")
        .attr("transform", "rotate(-90)")
        .text("Periods");

    var svgResiduals = d3.select("#coolStuff").append("svg")
        .attr("width", width)
        .attr("height", plotHeight)
        .append("g")
        .attr("class", "orbitSVG");

    var svgResYScale = d3.scale.linear()
        .domain([minDelay, maxDelay])
        .range([plotHeight - padding, padding]); 

    var svgResXScale = d3.scale.linear()
        .domain([0, tMax])
        .range([padding, width - padding]); 


    var svgResXAxis = d3.svg.axis()
        .scale(svgResXScale)
        .orient("bottom") // orientation of labels
        .ticks(5);  //Set rough # of ticks. D3 will take this as suggestion and decide based on "nice", whole intervals

    svgResiduals.append("g")
        .attr("class", "xaxis")  //Assign "axis" class
        .attr("transform", "translate(0," + (plotHeight - padding) + ")") // To align it to the bottom
        .call(svgResXAxis);

    var svgResYAxis = d3.svg.axis()
        .scale(svgResYScale)
        .orient("left") // orientation of labels
        .ticks(5);  //Set rough # of ticks. D3 will take this as suggestion and decide based on "nice", whole intervals

    svgResiduals.append("g")
        .attr("class", "yaxis")
        .attr("transform", "translate(" + padding + ",0)")
        .call(svgResYAxis);

    svgResiduals.append("text")
        .attr("class", "y label")
        .attr("text-anchor", "end")
        .attr("y", 6)
        .attr("dy", ".75em")
        .attr("transform", "rotate(-90)")
        .text("TOA delay");

    /************************
    * Define elliptical path - this and the clock ticks mechanism were heavily copied from
    * http://bl.ocks.org/cloudshapes/5662234
    *************************/

    var circleData = [];

    // Elliptical motion - pulsar
    t_circle = d3.map();
    t_circle.set("id", 2);
    t_circle.set("cr", 5);
    t_circle.set("rotrx", semimajAxis);
    t_circle.set("rotry", semiminAxis);
    t_circle.set("rtype", "ellipse");
    t_circle.set("offsetX", focusOffX);
    t_circle.set("offsetY", focusOffY);
    t_circle.set("timelimit", 2000);
    t_circle.set("starttime", undefined);
    t_circle.set("elapsed", 0);
    t_circle.set("x", 100);
    t_circle.set("y", 100);
    t_circle.set("scale", 1);
    t_circle.set("move", false);
    t_circle.set("angle", omegaPSkyRAD);
    t_circle.set("orbit_period", orbitPeriod);
    circleData.push(t_circle);


    var nextPulse = 0;

    var pulsations = [1];

    var ellipse = svg.append("ellipse")
    .attr("stroke-width", 2)
    .attr("rx", semimajAxis)
    .attr("ry", semiminAxis)
    .attr("fill", "none")
    .attr("stroke", "#ffffff")
    .attr("transform", function(d) {
      return "translate(" +(width / 2)  + ',' + (height / 2)  + ") rotate(" + omegaPSky +  ")" 
    });

    var majorAxis = svg.append("line")
    .attr("stroke-width", 2)
    .attr("x1", phi0_x)
    .attr("y1", phi0_y)
    .attr("x2", phi180_x)
    .attr("y2", phi180_y)
    .attr("stroke", "#ffffff");


    var T0Text = svg.append("text")
    .attr("class", "phaseLabel")
    .attr("text-anchor", "center")
    .attr("x", phi0_x)
    .attr("y", phi0_y)
    .text("T0");


    var minorAxis = svg.append("line")
    .attr("stroke-width", 2)
    .attr("x1", phi90_x)
    .attr("y1", phi90_y)
    .attr("x2", phi270_x)
    .attr("y2", phi270_y)
    .attr("stroke", "#ffffff");


    var T90Text = svg.append("text")
    .attr("class", "phaseLabel")
    .attr("text-anchor", "center")
    .attr("x", phi90_x)
    .attr("y", phi90_y)
    .text("T90");


    var xFocus = svg.append("line")
    .attr("stroke-width", 2)
    .attr("x1", (width / 2) + focusOffX - distANode_x)
    .attr("y1", (height / 2) + focusOffY)
    .attr("x2", (width / 2) + focusOffX + distDNode_x)
    .attr("y2", (height / 2) + focusOffY)
    .attr("stroke", "#bb3333");


    if (omegaPSky != 180){
      var ascendingNodeText = svg.append("text")
      .attr("class", "phaseLabel")
      .attr("text-anchor", "center")
      .attr("x", (width / 2) + focusOffX - distANode_x)
      .attr("y", (height / 2) + focusOffY)
      .text("Tasc");

      var decendingNodeText = svg.append("text")
      .attr("class", "phaseLabel")
      .attr("text-anchor", "center")
      .attr("x", (width / 2) + focusOffX + distDNode_x)
      .attr("y", (height / 2) + focusOffY)
      .text("Tdesc");
    }



    var focus = svg.append("rect")
    .attr("stroke-width", 2)
    .attr("width", 2)
    .attr("height", 2)
    .attr("x", (width / 2) + focusOffX - 1)
    .attr("y", (height / 2) + focusOffY - 1)
    .attr("fill", "red")
    .attr("stroke", "black");


    // Create the circle elements, and move them into their start positions
    var circle = svg.selectAll("circle")
    .data(circleData, function(d) { return d.get('id');})
    .enter()
    .append("circle")
    .attr("r", function(d)  { return d.get('cr'); })
    .attr('d', function(d) {
      d.set('x', phi0_x);
      d.set('y', phi0_y);
      return d;
    })
    .attr("transform", function(d) {return "translate(" + d.get('x') + "," + d.get('y') + "), scale(" + d.get('scale') + ")";});


    // timer_ret_val: could be used to stop the timer, but not actually used in this code really. 
    var timer_ret_val = false;

    // Keeps a record of the elapsed time since the timer began.
    var timer_elapsed = 0;

    // Kick off the timer, and the action begins: 
    d3.timer(tickFn);

    var lastpos = [0,0];
    var lastt = 0;
    var delays = [];
    var periods = [];

    var crossNode = false;
    var crossT90 = false;

    var ySave = "none";

    var firstLoop = true;

    function tickFn(_elapsed) {
      timer_elapsed = _elapsed;

      var t_circleData = circleData[0];


      if (t_circleData.get('move') == true) {
        if (t_circleData.get('starttime') == undefined){
          t_circleData.set('starttime', _elapsed);
        }

        // Calc elapsed time.
        var t_elapsed = _elapsed - t_circleData.get('starttime');

        // Keep a record.
        t_circleData.set('elapsed', t_elapsed);

        // Calculate how far through the desired time for one iteration.
        var t = t_elapsed / t_circleData.get('timelimit');

        var rotation_radius_x = t_circleData.get('rotrx');
        var rotation_radius_y = t_circleData.get('rotry');

        var t_offsetX = t_circleData.get('offsetX');
        var t_offsetY = t_circleData.get('offsetY');
        
        // Mean anomaly
        var t_angle = 1 / t_circleData.get('orbit_period') * t;
        t_angle -= Math.floor(t_angle);
        t_angle *= TWOPI;

        // Eccentricity
        var ecc = Math.sqrt(1 - Math.pow(t_circleData.get('rotry') / t_circleData.get('rotrx'), 2));

        // Eccentric anomaly
        var E = EccAnom(ecc, t_angle / TWOPI * 360.0, 5);

        // Position
        var pos = position(t_circleData.get('rotrx'), ecc, E / 360. * TWOPI);

        var t_x = pos[0];
        var t_y = pos[1];

        var dpos = Math.sqrt(Math.pow(lastpos[0] - pos[0], 2) + Math.pow(lastpos[1] - pos[1], 2))

        // rotate figure according to angle
        var dist = Math.sqrt(t_x * t_x + t_y * t_y);
        var old_angle = Math.acos(t_x / dist) * sign(t_y)
        var new_angle = t_circleData.get('angle') + old_angle;

        t_x = dist * Math.cos(new_angle);
        t_y = dist * Math.sin(new_angle);

        if (sign(t_y) != sign(ySave) && ySave != "none"){crossNode = true;};

        t_circleData.set('x', (width/2) + t_offsetX+ t_x);
        t_circleData.set('y', (height/2) + t_offsetY + t_y);

        var radialVelocity =  (t_y - lastpos[1]) / (t_elapsed - lastt);
        /* Finally, calculate delays and periods */
        var delay = -t_y;
        delays.push([t_elapsed, delay]);
        var obsPeriod = pulsePeriod * (1 + 0.001 * radialVelocity);
        periods.push([t_elapsed, obsPeriod]);

        if (delay > maxDelay && (!firstLoop)){maxDelay = delay};
        if (delay < minDelay && (!firstLoop)){minDelay = delay};
        if (obsPeriod > maxPeriod && (!firstLoop)){maxPeriod = obsPeriod};
        if (obsPeriod < minPeriod && (!firstLoop)){minPeriod = obsPeriod};

        if (t_elapsed > tMax && (!firstLoop)) {tMax=t_elapsed;};
        if (tMax - tMin > tSpan && (!firstLoop)) {tMin = tMax - tSpan;};

        ySave = t_y;


        lastpos = [t_x, t_y];
        lastt = t_elapsed;
        firstLoop = false;

      }

      // Actually move the circles and the text.
      var t_circle = svg.selectAll("circle");
      t_circle
      .attr("transform", function(d) {return "translate(" + d.get('x') + "," + d.get('y') + "), scale(" + d.get('scale') + ")";});

      svgResXScale.domain([0, tMax]);
      svgResYScale.domain([minDelay, maxDelay]);

      svgPerXScale.domain([0, tMax]);
      svgPerYScale.domain([minPeriod, maxPeriod]);

      svgResiduals.select(".xaxis")
                        .transition()  // https://github.com/mbostock/d3/wiki/Transitions#wiki-d3_ease
                        .call(svgResXAxis);   

      svgPeriods.select(".xaxis")
                        .transition()  // https://github.com/mbostock/d3/wiki/Transitions#wiki-d3_ease
                        .call(svgPerXAxis);   

      var delay_symb = svgResiduals.selectAll("circle").data(delays)
      
      var circleClass = crossNode == true ? "node" : "plot";
      var circleRadius = crossNode == true ? 4 : 1;

      crossNode = false;

      delay_symb.exit().remove();

      delay_symb.transition()
      .attr("cx", function(d){return svgResXScale(d[0]);})
      .attr("cy", function(d){return svgResYScale(d[1]);})

      delay_symb.enter()
      .append("svg:circle")
      .attr("r", circleRadius)
      .attr("cx", function(d){return svgResXScale(d[0]);})
      .attr("cy", function(d){return svgResYScale(d[1]);})
      .attr("class", circleClass);

      var per_symb = svgPeriods.selectAll("circle").data(periods)
      
      per_symb.exit().remove();

      per_symb.transition()
      .attr("cx", function(d){return svgPerXScale(d[0]);})
      .attr("cy", function(d){return svgPerYScale(d[1]);})

      per_symb.enter()
      .append("svg:circle")
      .attr("r", circleRadius)
      .attr("cx", function(d){return svgPerXScale(d[0]);})
      .attr("cy", function(d){return svgPerYScale(d[1]);})
      .attr("class", circleClass);

      return timer_ret_val;
    }


    // These two buttons don't stop the timer, 
    // just set flags to indicate that these two elements should stop/start:
    var startbtn=d3.select("#startbtn");
    startbtn.on("click", function() {
      for (var i = 0; i<circleData.length;i++)  {
        var t_circleData = circleData[i];
        t_circleData.set('move', true);
        t_circleData.set('starttime', timer_elapsed - t_circleData.get('elapsed'));
      }});


    var stopbtn=d3.select("#stopbtn");
    stopbtn.on("click", function()  {
      for (var i = 0; i<circleData.length;i++)  {
        var t_circleData = circleData[i];
        t_circleData.set('move', false);
      }});
} // end of doAll

var refrbtn=d3.select("#refrbtn");
refrbtn.on("click", function()  {
  d3.selectAll("svg")
       .remove();
       
  doAll();

});

