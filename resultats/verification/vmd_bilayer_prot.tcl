# Verification des feuillets de bilayer_prot dans VMD
# Bleu / orange : feuillets de REFERENCE (leaflets.txt) ; rouge : lipides que V4 place AUTREMENT.
# Lancer depuis le dossier du projet : vmd -e resultats/verification/vmd_bilayer_prot.tcl
mol new {Projet3/bilayer_prot/bilayer_prot.gro} waitfor all
mol delrep 0 top
mol representation Points 3
mol color ColorID 0
mol selection {residue 160 161 162 163 164 165 166 167 168 169 170 171 172 173 189 190 191 192 193 194 195 196 197 198 199 200 201 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 245 246 247 248 249 250 251 252 253 254 255 256 257 258}
mol addrep top
mol representation Points 3
mol color ColorID 3
mol selection {residue 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 232 233 234 235 236 237 238 239 240 241 242 243 244 259 260 261 262 263 264 265 266 267 268 269 270 271}
mol addrep top
mol representation VDW 1.4 12
mol color ColorID 1
mol selection {none}
mol addrep top
display projection Orthographic
axes location off
catch {package require pbctools; pbc box}
puts "bilayer_prot : 0 lipides mal places par V4 (en rouge)"
