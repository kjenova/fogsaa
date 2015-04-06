immutable AlignmentScoring
	match::Int
	mismatch::Int
	gap_start::Int
	gap::Int
end

immutable FOGSAABranch
	m::Int
	x::Int
	y::Int
end

function FOGSAA(X, Y, s::AlignmentScoring)
	const XY = 0x1
	const xY = 0x2
	const Xy = 0x3

	g0 = endof(X) - endof(Y)
	gM = g0 == 0? 0 : s.gap_start

	if g0 < 0
		gM += -g0 * s.gap		
		m0 = endof(X)
	else
		gM += g0 * s.gap
		m0 = endof(Y)
	end

	gm = gM + m0 * s.mismatch
	gM += m0 * s.match

	alignments = fill(0x0, endof(X) + 1, endof(Y))
	scores = fill(gm, endof(X) + 1, endof(Y) + 1)
	scores[1, 1] = 0

	Mi = mi = gM - gm + 1
	queue = Vector{Vector{FOGSAABranch}}(Mi)
	queue[Mi] = [FOGSAABranch(1, 1, 1)]

	ox = endof(X) + 1
	oy = endof(Y) + 1

	function enqueue(M, m, x, y, alignment)
		if alignment & 1 != 0
			x = -x
		end

		if alignment & 2 != 0
			y = -y
		end

		branch = FOGSAABranch(m, x, y)

		i = M - gm + 1

		if isdefined(queue, i)
			j = 1

			while j <= endof(queue[i]) && queue[i][j].m <= m
				j += 1
			end

			insert!(queue[i], j, branch)
		else
			queue[i] = [branch]
		end
	end

	function dequeue()
		x = queue[Mi][end].x
		y = queue[Mi][end].y

		pop!(queue[Mi])

		alignment = 0x0

		if x < 0
			x = -x
			alignment |= 1
		end

		if y < 0
			y = -y
			alignment |= 2
		end

		while !isdefined(queue, Mi) || endof(queue[Mi]) == 0
			if Mi == mi
				Mi = mi = 0
				break
			end

			Mi -= 1
		end

		x, y, alignment
	end

	while Mi - 1 + gm > scores[ox, oy]
		x, y, alignment = dequeue()

		while x <= endof(X) && y <= endof(Y)
			sorted = 0x0
			score = [(scores[x, y] +
				  (X[x] == Y[y]? s.match : s.mismatch)),
				 (scores[x, y] + s.gap +
				  (alignment == xY? 0 : s.gap_start)),
				 (scores[x, y] + s.gap +
				  (alignment == Xy? 0 : s.gap_start))]
			Mm = Vector{Int}(6)

			if score[XY] > scores[x + 1, y + 1]
				scores[x + 1, y + 1] = score[XY]
			let
				rx = endof(X) - x
				ry = endof(Y) - y
				g = rx - ry
				m = ry

				Mm[2*XY-1] = score[XY]

				if g != 0
					Mm[2*XY-1] += s.gap_start

					if g < 0
						g = -g
						m = rx
					end
				end

				Mm[2*XY-1] += g * s.gap
				Mm[2*XY] = Mm[2*XY-1] + m * s.mismatch
				Mm[2*XY-1] += m * s.match

				if Mm[2*XY-1] > scores[ox, oy]
					sorted = XY
				end
			end
			end

			if score[xY] > scores[x, y + 1]
				scores[x, y + 1] = score[xY]
			let
				rx = endof(X) + 1 - x
				ry = endof(Y) - y
				g = ry - rx
				m = rx

				Mm[2*xY-1] = score[xY]

				if g < 0
					g = -g
					m = ry
					Mm[2*xY-1] += s.gap_start
				end

				Mm[2*xY-1] += g * s.gap
				Mm[2*xY] = Mm[2*xY-1] + m * s.mismatch
				Mm[2*xY-1] += m * s.match

				if Mm[2*xY-1] > scores[ox, oy]
					if sorted == 0x0 ||
					  (Mm[2*xY-1] >= Mm[2*XY-1] &&
					   (Mm[2*xY-1] > Mm[2*XY-1] ||
					    Mm[2*xY] > Mm[2*XY]))
						sorted <<= 2
						sorted |= xY
					else
						sorted |= xY << 2
					end
				end
			end
			end

			if score[Xy] > scores[x + 1, y]
				scores[x + 1, y] = score[Xy]
			let
				rx = endof(X) - x
				ry = endof(Y) + 1 - y
				g = rx - ry
				m = ry

				Mm[2*Xy-1] = score[Xy]

				if g < 0
					g = -g
					m = rx
					Mm[2*Xy-1] += s.gap_start
				end

				Mm[2*Xy-1] += g * s.gap
				Mm[2*Xy] = Mm[2*Xy-1] + m * s.mismatch
				Mm[2*Xy-1] += m * s.match

				if Mm[2*Xy-1] > scores[ox, oy]
					t = (sorted & 0b1100) >> 2

					if t == 0 ||
					  (Mm[2*Xy-1] >= Mm[2*t-1] &&
					   (Mm[2*Xy-1] > Mm[2*t-1] ||
					    Mm[2*Xy] > Mm[2*t]))
						u = sorted & 0b11

						if u == 0 ||
						  (Mm[2*Xy-1] >= Mm[2*u-1] &&
						   (Mm[2*Xy-1] > Mm[2*u-1] ||
						    Mm[2*Xy] > Mm[2*u]))
							sorted <<= 2
							sorted |= Xy
						else
							sorted &= 0b11
							t = t << 4 | Xy << 2
							sorted |= t
						end
					else
						sorted |= Xy << 4
					end
				end
			end
			end

			if sorted == 0x0
				@goto prune
			end

			alignment = sorted & 0b11
			x = alignment & (XY | Xy) == 0? x : x + 1
			y = alignment & (XY | xY) == 0? y : y + 1

			a = (sorted & 0b1100) >> 2

			if a != 0 && Mm[2*a-1] > Mm[2*alignment]
				x1 = a & (XY | Xy) == 0? x : x + 1
				y1 = a & (XY | xY) == 0? y : y + 1

				enqueue(Mm[2*a-1], Mm[2*a], x1, y1, a)

				b = (sorted & 0b110000) >> 4

				if b != 0 && Mm[2*b-1] > Mm[2*a]
					x2 = b & (XY | Xy) == 0? x : x + 1
					y2 = b & (XY | xY) == 0? y : y + 1

					enqueue(Mm[2*b-1], Mm[2*b], x2, y2, b)
				end
			end
		end

		ox = x
		oy = y
@label prune
	end
end
