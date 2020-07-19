home=/Users/zhaomingyin/gitlocal/optkit/

#dot -Tpng vis_adj.dot > vis_adj.png
#dot -Tpng vis_tmp.dot > vis_tmp.png
#dot -Tpng vis_encode.dot > vis_encode.png
#dot -Tpng vis_a.dot > vis_a.png


mv ${home}vis_csr.dot  ${home}vis/vis_csr.dot
mv ${home}vis_adj.dot  ${home}vis/vis_adj.dot
mv ${home}vis_a.dot  ${home}vis/vis_a.dot
mv ${home}vis_encode.dot  ${home}vis/vis_encode.dot

dot -Tsvg ${home}vis/vis_csr.dot > ${home}vis/vis_csr.svg
dot -Tsvg ${home}vis/vis_adj.dot > ${home}vis/vis_adj.svg
dot -Tsvg ${home}vis/vis_a.dot > ${home}vis/vis_a.svg
dot -Tsvg ${home}vis/vis_encode.dot > ${home}vis/vis_encode.svg
