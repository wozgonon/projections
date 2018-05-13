

DIR=../target

MERCATOR=../target/debug/examples/geojson_to_web_mercator
DRAW=../examples/draw_to_js.py

all: scotland bouldercolerado

scotland: ${DIR}/scotland.html

bouldercolerado: ${DIR}/bouldercolerado.html

${DIR}/scotland.html: ${DIR}/wpc.json
	${MERCATOR} $< | ${DRAW} > $@

${DIR}/bouldercolerado.html: ${DIR}/bouldercolerado.geo.json
	${MERCATOR} $< | ${DRAW} > $@

world: ${DIR}/world.html

${DIR}/world.html: ../data/custom.geo.json
	${MERCATOR} $< | ${DRAW} > $@

${DIR}/wpc.json:
	curl --silent -o $@ https://raw.githubusercontent.com/martinjc/UK-GeoJSON/master/json/electoral/sco/wpc.json

${DIR}/bouldercolerado.geo.json:
	curl --silent -o $@ https://www-static.bouldercolorado.gov/docs/opendata/Subdivisions.GeoJSON


#https://geojson-maps.ash.ms/

clean:
	rm -f ${DIR}/*.html ${DIR}/*.json
