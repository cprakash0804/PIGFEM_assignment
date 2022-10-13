/*
 * head.h
 *
 *  Created on: Apr 16, 2019
 *      Author: reza
 */

#ifndef HEAD_H_
#define HEAD_H_

bool penetration_store(bool a)
{
	static bool p = 0;
	if (a==0)
	{
		return p;
	}else{
		if(p==0)
		{
			p=1;
		}else{
			p=0;
		}
		return p;
	}

}

#endif /* HEAD_H_ */
